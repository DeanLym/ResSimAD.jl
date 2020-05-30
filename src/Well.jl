module Well

using DataFrames: DataFrame

using ..Global: α
using ..AutoDiff:Tensor, zeros_tensor
using ..Grid:AbstractStructGrid, get_grid_index
using ..State:AbstractState

abstract type AbstractWell end

@enum CtrlMode SHUT CBHP CORAT CWRAT CGRAT CLRAT
@enum Limit MAX_BHP MIN_BHP MAX_ORAT MAX_WRAT MAX_GRAT MAX_LRAT
@enum WellType PRODUCER INJECTOR

const get_ctrl_mode = Dict{String,CtrlMode}(
    "bhp" => CBHP,
    "orat" => CORAT,
    "wrat" => CWRAT,
    "grat" => CGRAT,
    "lrat" => CLRAT,
)

mutable struct StandardWell{T} <: AbstractWell
    name::String
    r::Float64 # Wellbore radius
    ind::Vector{Int} # Perforation gridblock indices
    wi::Vector{Float64} # Well indices

    mode::CtrlMode
    target::Float64

    limits::Dict{Limit, Float64}

    qo::Tensor # Oil Rate
    qw::Tensor # Water Rate
    pw::Tensor # Wellbore pressure

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
    numvar::Int,
) where {T}
    p = StandardWell{T}(name)
    p.ind = Vector{Int}()
    append!(p.ind, perforation)

    p.r = radius

    nperf = length(p.ind) # Number of perforations
    p.wi = Vector{Float64}(undef, nperf)

    p.qo = zeros_tensor(nperf, numvar)
    p.qw = zeros_tensor(nperf, numvar)
    p.pw = zeros_tensor(nperf, numvar)

    return p
end

function init_results_df()
    Vec = Vector{Float64}
    DataFrame(Time=Vec(), ORAT=Vec(), WRAT=Vec(), GRAT=Vec(), WBHP=Vec())
end

function compute_wi(well::AbstractWell, grid::AbstractStructGrid)::Nothing
    ind, r = well.ind, well.r
    kx, ky, dx, dy, dz = grid.kx, grid.ky, grid.dx, grid.dy, grid.dz
    for (i, n) in enumerate(ind)
        kx, ky, dx, dy, dz = kx[n], ky[n], dx[n], dy[n], dz[n]
        r0 = 0.28*√(√(ky/kx)*dx^2 + √(kx/ky)*dy^2) / ((ky/kx)^0.25 + (kx/ky)^0.25)
        well.wi[i] = 2*α*pi*√(kx*ky)*dz / log(r0 / r)
    end
    return nothing
end





end
