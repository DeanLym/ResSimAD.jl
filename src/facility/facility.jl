module Facility

using DataFrames: DataFrame
using Statistics: mean

using Memento
const LOGGER = getlogger(@__MODULE__)
__init__() = Memento.register(LOGGER)


using ..Global: α
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
    p.r = radius
    nperf = length(p.ind) # Number of perforations
    p.wi = Vector{Float64}(undef, nperf)

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
    nperf = length(p.ind) # Number of perforations
    p.wi = Vector{Float64}(undef, nperf)

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
    bhp = mean(value(well.bhp))
    push!(well.results, [t, qo, qw, qg, qo+qw, bhp])
end

function compute_wi(well::StandardWell, grid::AbstractStructGrid, rock::AbstractRock)::Nothing
    ind, r = well.ind, well.r
    kx, ky, dx, dy, dz = rock.kx, rock.ky, grid.dx, grid.dy, grid.dz
    for (i, n) in enumerate(ind)
        kx, ky, dx, dy, dz = kx[n], ky[n], dx[n], dy[n], dz[n]
        r0 = 0.28*√(√(ky/kx)*dx^2 + √(kx/ky)*dy^2) / ((ky/kx)^0.25 + (kx/ky)^0.25)
        well.wi[i] = 2*α*pi*√(kx*ky)*dz / log(r0 / r)
    end
    return nothing
end

function compute_qo(well::StandardWell{PRODUCER}, fluid::AbstractFluid)::ADVector
    mode, target, ind = well.mode, well.target, well.ind
    o, w = fluid.phases.o, fluid.phases.w
    if mode == CBHP
        λo, po = o.λ, o.p
        if any(po[ind] .< target)
            @. well.qo = 0.0 * po[ind]
        else
            @. well.qo = well.wi * λo[ind] * (po[ind] - target)
        end
    elseif mode == CORAT
        @. well.qo = target + 0.0 * o.p[ind]
    elseif mode == CLRAT
        λo, λw = o.λ, w.λ
        @. well.qo = λo[ind] / (λo[ind] + λw[ind]) * target
    elseif mode == SHUT
        @. well.qo = 0.0 * o.p[ind]
    end
    return well.qo
end

function compute_qw(well::StandardWell{PRODUCER}, fluid::AbstractFluid)::ADVector
    mode, target, ind = well.mode, well.target, well.ind
    o, w = fluid.phases.o, fluid.phases.w
    if mode == CBHP
        λw, pw = w.λ, w.p
        if any(pw[ind] .< target)
            @. well.qw = 0.0*pw[ind]
        else
            @. well.qw = well.wi * λw[ind] * (pw[ind] - target)
        end
    elseif mode == CORAT
        λo, λw = o.λ, w.λ
        @. well.qw = (λw[ind] / λo[ind]) * target
    elseif mode == CLRAT
        λo, λw = o.λ, w.λ
        @. well.qw = λw[ind] / (λo[ind] + λw[ind]) * target
    elseif mode == SHUT
        @. well.qw = 0.0 * w.p[ind]
    else
        throw(ErrorException("Mode $mode not supported for producer"))
    end
    return well.qw
end


function compute_qo(well::StandardWell{INJECTOR}, fluid::AbstractFluid)
    @. well.qo = 0.0 * fluid.phases.o.p[well.ind]
    return well.qo
end

function compute_qw(well::StandardWell{INJECTOR}, fluid::AbstractFluid)::ADVector
    mode, target, ind = well.mode, well.target, well.ind
    o, w = fluid.phases.o, fluid.phases.w
    if mode == CBHP
        kro, krw, μo, μw, bw, pw = o.kr, w.kr, o.μ, w.μ, w.b, w.p
        λ = @. (kro[ind] / μo[ind] + krw[ind] / μw[ind]) / bw[ind]
        if any(pw[ind] .> target)
            # Avoid inverse flow
            @. well.qw = 0.0*pw[ind]
        else
            @. well.qw = well.wi * λ * (pw[ind] - target)
        end
    elseif mode == CWRAT
        @. well.qw = target + 0.0*w.p[ind]
    elseif mode == SHUT
        @. well.qw = 0.0*w.p[ind]
    else
        throw(ErrorException("Mode $mode not supported for injector"))
    end
    return well.qw
end


function compute_bhp(well::StandardWell{PRODUCER}, fluid::AbstractFluid)::ADVector
    mode, target, ind = well.mode, well.target, well.ind
    o, w = fluid.phases.o, fluid.phases.w
    if mode == CBHP
        @. well.bhp = target + 0.0*o.p[ind]
    elseif mode == CORAT
        λo, po = o.λ, o.p
        @. well.bhp = po[ind] - target / (well.wi * λo[ind])
    elseif mode == CLRAT
        λo, λw = o.λ, w.λ
        @. well.bhp = o.p[ind] - target / (well.wi * (λo[ind] + λw[ind]))
    elseif mode == SHUT
        @. well.bhp = o.p[ind]
    end
end

function compute_bhp(well::StandardWell{INJECTOR}, fluid::AbstractFluid)::ADVector
    mode, target, ind = well.mode, well.target, well.ind
    o, w = fluid.phases.o, fluid.phases.w

    if mode == CBHP
        @. well.bhp = target + 0.0*o.p[ind]
    elseif mode == CWRAT
        λo, λw, po = o.λ, w.λ, o.p
        @. well.bhp = po[ind] - target / (well.wi * (λo[ind] + λw[ind]))
    elseif mode == SHUT
        @. well.bhp = o.p[ind]
    end
end



function compute_ql(well::StandardWell, fluid::SPFluid)::ADVector
    mode, target, ind = well.mode, well.target, well.ind
    phase = fluid.phases[1]
    if mode == CBHP
        λ, p = phase.λ, phase.p
        @. well.ql = well.wi * λ[ind] * (p[ind] - target)
    elseif mode == CLRAT
        @. well.ql = target + 0.0 * phase.p[ind]
    elseif mode == SHUT
        @. well.ql = 0.0 * phase.p[ind]
    end
    return well.ql
end

include("check_limits.jl")

end
