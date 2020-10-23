module Facility

using DataFrames: DataFrame
using Statistics: mean

using Memento
const LOGGER = getlogger(@__MODULE__)
__init__() = Memento.register(LOGGER)


using ..Global: α, β, g_, gc
using ..AutoDiff:ADVector, advector, value
using ..Rock:AbstractRock
using ..Grid:AbstractStructGrid, get_grid_index
using ..Fluid:AbstractFluid, SPFluid

abstract type AbstractFacility end

@enum CtrlMode SHUT CBHP CORAT CWRAT CGRAT CLRAT
@enum Limit MAX_BHP MIN_BHP MAX_ORAT MAX_WRAT MAX_GRAT MAX_LRAT
@enum WellType PRODUCER INJECTOR

const get_ctrl_mode = Dict{String,CtrlMode}(
    "shut" => SHUT,
    "bhp" => CBHP,
    "orat" => CORAT,
    "wrat" => CWRAT,
    "grat" => CGRAT,
    "lrat" => CLRAT,
)

const get_limit = Dict{String, Limit}(
    "max_bhp" => MAX_BHP,
    "min_bhp" => MIN_BHP,
    "max_orat" => MAX_ORAT,
    "max_wrat" => MAX_WRAT,
    "max_grat" => MAX_GRAT,
    "max_lrat" => MAX_LRAT,
)

mutable struct StandardWell{T} <: AbstractFacility
    name::String
    r::Float64 # Wellbore radius
    ind::Vector{Int} # Perforation gridblock indices
    wi::Vector{Float64} # Well indices
    d::Vector{Float64} # Perforation depth
    ρl::Vector{Float64} # Liquid density

    mode::CtrlMode
    target::Float64

    limits::Dict{Limit, Float64}

    qo::ADVector # Oil Rate
    qw::ADVector # Water Rate
    ql::ADVector # Liquid Rate
    bhp::ADVector # Bottom hole pressure

    results::DataFrame

    function StandardWell{T}(name::String) where {T}
        @assert isa(T, WellType)
        x = new{T}()
        x.name = name
        x.limits = Dict{Limit, Float64}()
        x.results = init_results_df()
        x.mode = SHUT
        return x
    end
end

function StandardWell{T}(
    name::String,
    perforation::Vector{Int},
    radius::Float64,
    nv::Int,
) where {T}
    p = StandardWell{T}(name)
    p.ind = Vector{Int}()
    append!(p.ind, perforation)
    sort!(p.ind)
    p.r = radius
    nperf = length(p.ind) # Number of perforations
    p.wi = Vector{Float64}(undef, nperf)
    p.d = Vector{Float64}(undef, nperf)
    p.ρl = Vector{Float64}(undef, nperf)

    p.qo = advector(nperf, 1, nv)
    p.qw = advector(nperf, 1, nv)
    p.ql = advector(nperf, 1, nv)
    p.bhp = advector(nperf, 1, nv)

    return p
end

function StandardWell{T}(
    name::String,
    perforation::Vector{Int},
    nv::Int,
) where {T}
    p = StandardWell{T}(name)
    p.ind = Vector{Int}()
    append!(p.ind, perforation)
    sort!(p.ind)
    nperf = length(p.ind) # Number of perforations
    p.wi = Vector{Float64}(undef, nperf)
    p.d = Vector{Float64}(undef, nperf)
    p.ρl = Vector{Float64}(undef, nperf)

    p.qo = advector(nperf, 1, nv)
    p.qw = advector(nperf, 1, nv)
    p.ql = advector(nperf, 1, nv)
    p.bhp = advector(nperf, 1, nv)

    return p
end

function isproducer(::StandardWell{PRODUCER})
    return true
end

function isproducer(::StandardWell{INJECTOR})
    return false
end

function init_results_df()
    Vec = Vector{Float64}
    DataFrame(TIME=Vec(), ORAT=Vec(), WRAT=Vec(), GRAT=Vec(), LRAT=Vec(), WBHP=Vec())
end

function save_result(well::StandardWell, t::Float64)
    qo, qw, qg = sum(value(well.qo)), sum(value(well.qw)), 0.0
    bhp = value(well.bhp[1])
    push!(well.results, [t, qo, qw, qg, qo+qw, bhp])
end

function compute_wi(well::StandardWell, grid::AbstractStructGrid, rock::AbstractRock)::Nothing
    ind, r = well.ind, well.r
    kx_, ky_, dx_, dy_, dz_ = rock.kx, rock.ky, grid.dx, grid.dy, grid.dz
    for (i, n) in enumerate(ind)
        kx, ky, dx, dy, dz = kx_[n], ky_[n], dx_[n], dy_[n], dz_[n]
        r0 = 0.28*√(√(ky/kx)*dx^2 + √(kx/ky)*dy^2) / ((ky/kx)^0.25 + (kx/ky)^0.25)
        well.wi[i] = 2*α*pi*√(kx*ky)*dz / log(r0 / r)
    end
    return nothing
end

function compute_well_state(well::StandardWell{PRODUCER}, fluid::AbstractFluid)
    mode, target, ind, d, wi, ρl = well.mode, well.target, well.ind, well.d, well.wi, well.ρl
    o, w = fluid.phases.o, fluid.phases.w
    λo, po = o.λ, o.p
    λw, pw = w.λ, w.p
    ρo_, ρw_ = value(o.ρ[ind]), value(w.ρ[ind])
    λo_, λw_ = value(λo[ind]), value(λw[ind])
    # Compute well bore density   
    @. ρl = (ρo_ * λo_ + ρw_ * λw_)  / (λo_ + λw_ + 1.e-3)
    if mode == CBHP
        for (i, idx) in enumerate(ind)
            if po[idx] < target
                well.qo[i] = 0.0*po[idx]
            else
                well.qo[i] = wi[i] * λo[idx] * (po[idx] - target)
            end
            if pw[idx] < target
                well.qw[i] = 0.0*pw[idx]
            else
                well.qw[i] = wi[i] * λw[idx] * (pw[idx] - target)
            end
            well.bhp[i] = target + 0.0*o.p[idx]
            if i < length(ind)
                target += β / gc * g_ * ρl[i] * (d[i+1] - d[i])
            end
        end
    elseif mode == CORAT
        @. well.qo = target + 0.0 * po[ind]
        @. well.qw = (λw[ind] / λo[ind]) * target
        @. well.bhp = po[ind] - target / (wi * λo[ind])
        # rhs, lhs, Δpw = 0.0, 0.0, 0.0
        # po_ = value(po[ind])
        # for (i, idx) in enumerate(ind)
        #     lhs += wi[i] * λo_[i]
        #     rhs += wi[i] * λo_[i] * (po_[i] - Δpw)
        #     if i < length(ind)
        #         Δpw += β / gc * g_ * ρl[i] * (d[i+1] - d[i])
        #     end
        # end
        # Δpw = 0.0
        # for (i, idx) in enumerate(ind)
        #     well.bhp[i] = (rhs - target) / lhs + Δpw + 0.0 * po[idx]
        #     if i < length(ind)
        #         Δpw += β / gc * g_ * ρl[i] * (d[i+1] - d[i])
        #     end
        # end
        # @. well.qo = wi * λo[ind] * (po[ind] - well.bhp)
        # @. well.qw = wi * λw[ind] * (pw[ind] - well.bhp)
    elseif mode == CLRAT
        @. well.qo = λo[ind] / (λo[ind] + λw[ind]) * target
        @. well.qw = λw[ind] / (λo[ind] + λw[ind]) * target
        @. well.bhp = po[ind] - target / (wi * (λo[ind] + λw[ind]))
        # rhs, lhs, Δpw = 0.0, 0.0, 0.0
        # po_ = value(po[ind])
        # for (i, idx) in enumerate(ind)
        #     lhs += wi[i] * (λo_[i] + λw_[i])
        #     rhs += wi[i] * (λo_[i] + λw_[i]) * (po_[i] - Δpw)
        #     # sum{i}(wi[i] * (λo_[i] + λw_[i]) * (po_[i] - bhp[1] - ρl[i]gh[i])) = target
        #     # sum{i}(wi[i] * (λo_[i] + λw_[i]) * bhp[1]) = sum{i}(wi[i] * (λo_[i] + λw_[i]) * (po_[i] - ρl[i]gh[i])) - target
        #     # lhs * bhp[1] = rhs - target
        #     # bhp[1] = (rhs - target) / lhs
        #     if i < length(ind)
        #         Δpw += β / gc * g_ * ρl[i] * (d[i+1] - d[i])
        #     end
        # end
        # Δpw = 0.0
        # for (i, idx) in enumerate(ind)
        #     well.bhp[i] = (rhs - target) / lhs + Δpw + 0.0 * po[idx]
        #     if i < length(ind)
        #         Δpw += β / gc * g_ * ρl[i] * (d[i+1] - d[i])
        #     end
        # end
        # @. well.qo = wi * λo[ind] * (po[ind] - well.bhp)
        # @. well.qw = wi * λw[ind] * (pw[ind] - well.bhp)
    elseif mode == SHUT
        @. well.qo = 0.0 * po[ind]
        @. well.qw = 0.0 * pw[ind]
        @. well.bhp = po[ind]
    end
end

function compute_well_state(well::StandardWell{INJECTOR}, fluid::AbstractFluid)
    mode, target, ind, d, wi, ρl = well.mode, well.target, well.ind, well.d, well.wi, well.ρl
    o, w = fluid.phases.o, fluid.phases.w
    kro, krw, μo, μw, bw, pw, po = o.kr, w.kr, o.μ, w.μ, w.b, w.p, o.p
    @. well.qo = 0.0 * po[ind]
    @. ρl = value(w.ρ[ind])
    λ = @. (kro[ind] / μo[ind] + krw[ind] / μw[ind]) / bw[ind]
    if mode == CBHP
        for (i, idx) in enumerate(ind)
            if pw[idx] > target
                well.qw[i] = 0.0*pw[idx]
            else
                well.qw[i] = wi[i] * λ[i] * (pw[idx] - target)
            end
            well.bhp[i] = target + 0.0*w.p[idx]
            if i < length(ind)
                target += β / gc * g_ * ρl[i] * (d[i+1] - d[i])
            end
        end
    elseif mode == CWRAT
        @. well.qw = target + 0.0*pw[ind]
        @. well.bhp = po[ind] - target / (wi * λ)
        # λ_ = value(λ)
        # rhs, lhs, Δpw = 0.0, 0.0, 0.0
        # pw_ = value(pw[ind])
        # for (i, idx) in enumerate(ind)
        #     lhs += wi[i] * λ_[i]
        #     rhs += wi[i] * λ_[i] * (pw_[i] - Δpw)
        #     if i < length(ind)
        #         Δpw += β / gc * g_ * ρl[i] * (d[i+1] - d[i])
        #     end
        # end
        # Δpw = 0.0
        # for (i, idx) in enumerate(ind)
        #     well.bhp[i] = (rhs - target) / lhs + Δpw + 0.0 * pw[idx]
        #     if i < length(ind)
        #         Δpw += β / gc * g_ * ρl[i] * (d[i+1] - d[i])
        #     end
        # end
        # @. well.qw = wi * λ * (pw[ind] - well.bhp)
    elseif mode == SHUT
        @. well.qw = 0.0*pw[ind]
        @. well.bhp = po[ind]
    else
        throw(ErrorException("Mode $mode not supported for injector"))
    end
end

include("check_limits.jl")

end
