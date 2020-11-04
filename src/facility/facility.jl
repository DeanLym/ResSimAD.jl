module Facility

using DataFrames: DataFrame
using Statistics: mean

using Memento
const LOGGER = getlogger(@__MODULE__)
__init__() = Memento.register(LOGGER)


using ..Global: α, β, g_, gc
using ..AutoDiff:ADVector, advector, value, ADScaler, adscaler, grad
using ..Rock:AbstractRock
using ..Grid:AbstractStructGrid, get_grid_index
using ..Fluid:AbstractFluid

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

    row_num::Int # row number of well equation
    bhp::ADScaler # Bottom hole pressure

    qo::ADVector # Oil Rate
    qw::ADVector # Water Rate
    ql::ADVector # Liquid Rate

    ∂r::Array{Float64, 2} # Derivative to reservoir primary variables 
    ∂w::Float64 # Derivative to the well primary variable (bhp)  
    rw::Float64 # Residual of well equation

    bhp_::ADVector 

    
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
    nv::Int,
    row_num::Int,
) where {T}
    p = StandardWell{T}(name)
    p.ind = Vector{Int}()
    append!(p.ind, perforation)
    sort!(p.ind)
    nperf = length(p.ind) # Number of perforations
    p.wi = Vector{Float64}(undef, nperf)
    p.d = Vector{Float64}(undef, nperf)
    p.ρl = Vector{Float64}(undef, nperf)

    p.row_num = row_num

    p.bhp = adscaler(0.0, row_num, 1, nv)

    p.qo = advector(nperf, 2, nv)
    p.qw = advector(nperf, 2, nv)
    p.ql = advector(nperf, 2, nv)
    p.bhp_ = advector(nperf, 2, nv)

    p.∂r = zeros((nperf, nv))
    p.∂w = 0.0
    p.rw = 0.0

    return p
end

function StandardWell{T}(
    name::String,
    perforation::Vector{Int},
    radius::Float64,
    nv::Int,
    row_num::Int,
) where {T}
    p = StandardWell{T}(name, perforation, nv, row_num)
    p.r = radius
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
    bhp = value(well.bhp)
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
    bhp = well.bhp
    o, w = fluid.phases.o, fluid.phases.w
    λo, po = o.λ, o.p
    λw, pw = w.λ, w.p
    ρo_, ρw_ = value(o.ρ[ind]), value(w.ρ[ind])
    λo_, λw_ = value(λo[ind]), value(λw[ind])
    nv = size(well.∂r)[end]

    # Compute well rates
    if mode == SHUT
        # Well equation: rw = bhp - po
        @. well.qo = 0.0 * (po[ind] - bhp)
        @. well.qw = 0.0 * (pw[ind] - bhp)
    else
        # Compute well bore density   
        @. ρl = (ρo_ * λo_ + ρw_ * λw_)  / (λo_ + λw_ + 1.e-3)
        Δbhp = 0.0
        for (i, idx) in enumerate(ind)
            # if po[idx] < target
            #     well.qo[i] = 0.0*(po[idx] - bhp)
            # else
            #     well.qo[i] = wi[i] * λo[idx] * (po[idx] - (bhp + Δbhp))
            # end
            # if pw[idx] < target
            #     well.qw[i] = 0.0*(pw[idx] - bhp)
            # else
            #     well.qw[i] = wi[i] * λw[idx] * (pw[idx] - (bhp + Δbhp))
            # end
            well.qo[i] = wi[i] * λo[idx] * (po[idx] - (bhp + Δbhp))
            well.qw[i] = wi[i] * λw[idx] * (pw[idx] - (bhp + Δbhp))
            if i < length(ind)
                Δbhp += β / gc * g_ * ρl[i] * (d[i+1] - d[i])
            end
        end
    end

    if mode == CBHP
        # Well equation: rw = bhp - target
        well.rw = value(well.bhp) - target
        well.∂r .= 0.0
        well.∂w = 1.0
    elseif mode == CORAT
        # Well equation: rw = sum(well.qo) - target
        well.rw = sum(value.(well.qo)) - target
        well.∂w = sum(grad.(well.qo, 2, 1))
        for i = 1:nv
            well.∂r[:, i] .= grad.(well.qo, 1, i) # ∂r∂po
        end
    elseif mode == CLRAT
        # Well equation: rw = sum(well.qo + well.qw) - target
        well.rw = sum(value.(well.qo)) + sum(value.(well.qw)) - target
        well.∂w = sum(grad.(well.qo, 2, 1)) + sum(grad.(well.qw, 2, 1))
        for i = 1:nv
            well.∂r[:, i] .= grad.(well.qo, 1, i) .+ grad.(well.qw, 1, i) # ∂r∂po
        end
    elseif mode == SHUT
        # Well equation: rw = bhp - po[1]
        well.rw = value(well.bhp) - value(po[1])
        well.∂r .= 0.0
        well.∂w = 1.0
    end
end

function compute_well_state(well::StandardWell{INJECTOR}, fluid::AbstractFluid)
    mode, target, ind, d, wi, ρl = well.mode, well.target, well.ind, well.d, well.wi, well.ρl
    bhp = well.bhp
    o, w = fluid.phases.o, fluid.phases.w
    kro, krw, μo, μw, bw, pw, po = o.kr, w.kr, o.μ, w.μ, w.b, w.p, o.p
    @. well.qo = 0.0 * (po[ind] - bhp)
    @. ρl = value(w.ρ[ind])
    nv = size(well.∂r)[end]

    if mode == SHUT
        @. well.qw = 0.0* (pw[ind] - bhp)
    else
        λ = @. (kro[ind] / μo[ind] + krw[ind] / μw[ind]) / bw[ind]
        Δbhp = 0.0
        for (i, idx) in enumerate(ind)
            # if pw[idx] > target
            #     well.qw[i] = 0.0*(pw[idx] - bhp)
            # else
            #     # well.qw[i] = wi[i] * λw[idx] * (pw[idx] - (bhp + Δbhp))
            #     well.qw[i] = wi[i] * λ[i] * (pw[idx] - (bhp + Δbhp))
            # end
            well.qw[i] = wi[i] * λ[i] * (pw[idx] - (bhp + Δbhp))
            # well.bhp_[i] = target + 0.0*w.p[idx]
            if i < length(ind)
                Δbhp += β / gc * g_ * ρl[i] * (d[i+1] - d[i])
            end
        end
    end

    if mode == CBHP
        # Well equation: rw = bhp - target
        well.∂r .= 0.0
        well.∂w = 1.0
        well.rw = value(well.bhp) - target
    elseif mode == CWRAT
        # Well equation: rw = sum(qw) - target
        well.rw = sum(value.(well.qw)) - target
        well.∂w = sum(grad.(well.qw, 2, 1))
        for i = 1:nv
            well.∂r[:, i] .= grad.(well.qw, 1, i) # ∂r∂po
        end
    elseif mode == SHUT
        # Well equation: rw = bhp - pw[1]
        well.rw = value(well.bhp) - value(pw[1])
        well.∂r .= 0.0
        well.∂w = 1.0
    else
        throw(ErrorException("Mode $mode not supported for injector"))
    end
end

include("check_limits.jl")

end
