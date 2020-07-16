module Fluid

using DataFrames
using CSV
using Interpolations
using Interpolations:GriddedInterpolation
const Interp = GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}}

#! format: off
import Base
using ..Global: M, β, g_, gc
using ..AutoDiff: ADVector, advector, value, set_primary_variable
using ..Grid: ConnList

include("pvt.jl")
include("kr.jl")

#
struct Phase
    p::ADVector   # Cells Phase Pressure
    s::ADVector   # Cells Phase Saturation
    b::ADVector   # Cells Phase Formation Volume Factor
    μ::ADVector   # Cells Phase Viscosity
    kr::ADVector  # Cells Phase Relative Permeability
    λ::ADVector   # Cells Phase Mobility
    ρ::ADVector   # Cells Phase Density

    γ::ADVector   # Connections Phase Gravity Term
    ΔΨ::ADVector  # Connections Phase Potential Difference
    f::ADVector   # Connections Component Flux

    upstream_ind::Vector{Int} # Index of upstream cell for each connection
    ΔΨv::Vector{Float64} # Value of Δψ

    pn::Vector{Float64}  # Cells Phase Pressure t-1
    sn::Vector{Float64}  # Cells Phase Saturation t-1
    bn::Vector{Float64}  # Cells Phase Formation Volume Factor t-1

    ρs::Float64  # Phase density at standard condition
    pvt::AbstractPVT  # PVT function for formation and viscosity

    p_rec::Dict{Float64, Vector{Float64}}
    s_rec::Dict{Float64, Vector{Float64}}

    function Phase(nc::Int, nconn::Int, nv::Int, ρs::Float64, pvt::AbstractPVT)
        props = (:p, :s, :b, :μ, :kr, :λ, :ρ)
        tensor_cell = [advector(nc, 1, nv) for _ in props]
        props = (:γ, :ΔΨ, :f)
        tensor_conn = [advector(nconn, 2, nv) for _ in props]
        upstream_ind = zeros(Int64, nconn)
        ΔΨv = zeros(nconn)
        props = (:pn, :sn, :bn)
        vec_cell = [zeros(Float64, nc) for _ in props]
        props = (:p_rec, :s_rec)
        vec_rec = [Dict{Float64, Vector{Float64}}() for _ in props]
        return new(tensor_cell..., tensor_conn..., upstream_ind, ΔΨv, vec_cell..., ρs, pvt, vec_rec...)
    end
end


struct Component
    a::ADVector   # Cells Component Accumulation
    r::Vector{Float64}   # Cells Component Residual
    function Component(nc::Int, nv::Int)
        return new(advector(nc, 1, nv), zeros(nc))
    end
end

function compute_ρ(phase::Phase)::ADVector
    @. phase.ρ = phase.ρs / phase.b
end

function compute_λ(phase::Phase)::ADVector
    @. phase.λ = phase.kr / (phase.μ * phase.b)
end

function compute_b(phase::Phase)::ADVector where {T <: AbstractPVT}
    get_b(phase.pvt, phase.b, phase.p)
end

function compute_μ(phase::Phase)::ADVector where {T <: AbstractPVT}
    get_μ(phase.pvt, phase.μ, phase.p)
end

function compute_γ(phase::Phase, connlist::ConnList)::ADVector
    l, r = connlist.l, connlist.r
    @views ρl = phase.ρ[l]
    @views ρr = phase.ρ[r]
    @. phase.γ = β * g_ / gc * (ρl + ρr) / 2.0
end

function compute_ΔΨ(phase::Phase, connlist::ConnList)::ADVector
    l, r, Δd = connlist.l, connlist.r, connlist.Δd
    @views pl = phase.p[l]
    @views pr = phase.p[r]
    @. phase.ΔΨ = pl - pr - phase.γ * Δd
end

function get_upstream_ind(ΔΨv::Float64, l::Int, r::Int)
    if ΔΨv > 0.0
        return l
    else
        return r
    end
end

function compute_f(phase::Phase, connlist::ConnList)::ADVector
    l, r, trans = connlist.l, connlist.r, connlist.trans
    @. phase.ΔΨv = value(phase.ΔΨ)
    @. phase.upstream_ind = get_upstream_ind(phase.ΔΨv, l, r)
    @views λconn = phase.λ[phase.upstream_ind]
    @. phase.f = λconn * phase.ΔΨ * trans
end

function set_phase_tn(phase::Phase, p::Vector{Float64}, s::Vector{Float64})
    phase.pn .= p
    phase.sn .= s
    get_b(phase.pvt, phase.bn, p)
end

function update_phase_tn(phase::Phase)::Phase
    for (x, y) in zip((:pn, :sn, :bn), (:p, :s, :b))
        getfield(phase, x) .= value(getfield(phase, y))
    end
    return phase
end

function update_phase(phase::Phase, connlist::ConnList)::Phase
    # println("SGrad")
    # @time begin
    compute_b(phase)
    compute_μ(phase)
    compute_ρ(phase)
    compute_λ(phase)
    # end
    # println("DGrad")
    # @time begin
    compute_γ(phase, connlist)
    compute_ΔΨ(phase, connlist)
    compute_f(phase, connlist)
    # end

    return phase
end


## Define Abstract State
abstract type AbstractFluid end

function update_fluid_tn(fluid::AbstractFluid)
    for phase in fluid.phases
        update_phase_tn(phase)
    end
end

## OWFluid
struct OWFluid <: AbstractFluid
    nv::Int
    phases::NamedTuple{(:o, :w), Tuple{Phase, Phase}}
    components::NamedTuple{(:o, :w), Tuple{Component, Component}}
    krow::AbstractKROW
    function OWFluid(nc::Int, nconn::Int, ρo::Float64, ρw::Float64,
        pvto::AbstractPVT, pvtw::AbstractPVT, krow::AbstractKROW)
        nv = 2
        phases = (o = Phase(nc, nconn, nv, ρo, pvto), w=Phase(nc, nconn, nv, ρw, pvtw))
        components = (o = Component(nc, nv), w= Component(nc, nv))
        return new(nv, phases, components, krow)
    end
end

function compute_kr(fluid::OWFluid)::OWFluid
    phases = fluid.phases
    get_kro(fluid.krow, phases.o.kr, phases.o.s)
    get_krw(fluid.krow, phases.w.kr, phases.w.s)
    return fluid
end

function set_fluid_tn(fluid::OWFluid, po::Vector{Float64}, sw::Vector{Float64})
    phases = fluid.phases
    set_phase_tn(phases.o, po, 1 .- sw)
    set_phase_tn(phases.w, po, sw)
end

function update_primary_variable(fluid::OWFluid, po::Vector{Float64}, sw::Vector{Float64})
    phases = fluid.phases
    set_primary_variable(phases.o.p, po, 1)
    set_primary_variable(phases.w.s, sw, 2)
end

function reset_primary_variable(fluid::OWFluid)
    phases = fluid.phases
    update_primary_variable(fluid, phases.o.pn, phases.w.sn)
end

function update_phases(fluid::OWFluid, connlist::ConnList)
    @. fluid.phases.w.p = fluid.phases.o.p
    @. fluid.phases.o.s = 1 - fluid.phases.w.s
    # println("Compute Kr")
    compute_kr(fluid)
    update_phase(fluid.phases.o, connlist)
    update_phase(fluid.phases.w, connlist)
end

## Single Phase Fluid
struct SPFluid <: AbstractFluid
    nv::Int
    phases::NamedTuple
    components::NamedTuple
    function SPFluid(nc::Int, nconn::Int, ρ::Float64, pvt::AbstractPVT; name::String = "OIL")
        nv = 1
        if name == "OIL"
            phases = (o = Phase(nc, nconn, nv, ρ, pvt), )
            components = (o = Component(nc, nv), )
        elseif name == "WATER"
            phases = (w = Phase(nc, nconn, nv, ρ, pvt), )
            components = (w = Component(nc, nv), )
        elseif name == "GAS"
            phases = (g = Phase(nc, nconn, nv, ρ, pvt), )
            components = (g = Component(nc, nv), )
        end

        return new(nv, phases, components)
    end
end

function set_fluid_tn(fluid::SPFluid, p::Vector{Float64})
    phase = fluid.phases[1]
    @. phase.sn = 1.0
    set_phase_tn(phase, p, phase.sn)
end

function update_primary_variable(fluid::SPFluid, p::Vector{Float64})
    phase = fluid.phases[1]
    set_primary_variable(phase.p, p, 1)
    @. phase.kr = 1.0 + 0.0*phase.p
    @. phase.s = 1.0 + 0.0*phase.p
end

function reset_primary_variable(fluid::SPFluid)
    phase = fluid.phases[1]
    update_primary_variable(fluid, phase.pn)
end

function update_phases(fluid::SPFluid, connlist::ConnList)
    update_phase(fluid.phases[1], connlist)
end


end
