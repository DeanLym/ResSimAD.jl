module Facility

using DataFrames: DataFrame
using Statistics: mean

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
    pw::ADVector # Wellbore pressure

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
    p.pw = advector(nperf, 1, nv)

    return p
end

function init_results_df()
    Vec = Vector{Float64}
    DataFrame(TIME=Vec(), ORAT=Vec(), WRAT=Vec(), GRAT=Vec(), LRAT=Vec(), WBHP=Vec())
end

function save_result(well::StandardWell, t::Float64)
    qo, qw, qg, ql = sum(value(well.qo)), sum(value(well.qw)), 0.0, sum(value(well.ql))
    pw = mean(value(well.pw))
    push!(well.results, [t, qo, qw, qg, ql, pw])
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
        @. well.qo = well.wi * λo[ind] * (po[ind] - target)
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
        @. well.qw = well.wi * λw[ind] * (pw[ind] - target)
    elseif mode == CORAT
        λo, λw = o.λ, w.λ
        @. well.qw = (λw[ind] / λo[ind]) * target
    elseif mode == CLRAT
        λo, λw = o.λ, w.λ
        @. well.qw = λw[ind] / (λo[ind] + λw[ind]) * target
    elseif mode == SHUT
        @. well.qw = 0.0 * w.p[ind]
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
        λo, λw, pw = o.λ, w.λ, w.p
        @. well.qw = well.wi * (λo[ind] + λw[ind]) * (pw[ind] - target)
    elseif mode == CWRAT
        @. well.qw = target + 0.0*w.p[ind]
    elseif mode == SHUT
        @. well.qw = 0.0*w.p[ind]
    end
    return well.qw
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

end
