module Fluid

using DataFrames
using CSV
using Interpolations
using Interpolations:GriddedInterpolation
const Interp = GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}}

#! format: off
import Base
using ..Global: M, β, g_, gc
using ..AutoDiff: Tensor, param, zeros_tensor, data
using ..Grid: ConnList

include("pvt.jl")
include("kr.jl")

#
struct Phase
    p::Tensor   # Cells Phase Pressure
    s::Tensor   # Cells Phase Saturation
    b::Tensor   # Cells Phase Formation Volume Factor
    μ::Tensor   # Cells Phase Viscosity
    kr::Tensor  # Cells Phase Relative Permeability
    λ::Tensor   # Cells Phase Mobility
    ρ::Tensor   # Cells Phase Density

    γ::Tensor   # Connections Phase Gravity Term
    ΔΨ::Tensor  # Connections Phase Potential Difference
    f::Tensor   # Connections Component Flux

    pn::Vector{Float64}  # Cells Phase Pressure t-1
    sn::Vector{Float64}  # Cells Phase Saturation t-1
    bn::Vector{Float64}  # Cells Phase Formation Volume Factor t-1

    ρs::Float64  # Phase density at standard condition
    pvt::AbstractPVT  # PVT function for formation and viscosity

    p_rec::Dict{Float64, Vector{Float64}}
    s_rec::Dict{Float64, Vector{Float64}}

    function Phase(nc::Int, nconn::Int, nv::Int, ρs::Float64, pvt::AbstractPVT)
        props = (:p, :s, :b, :μ, :kr, :λ, :ρ)
        tensor_cell = [zeros_tensor(nc, nv) for _ in props]
        props = (:γ, :ΔΨ, :f)
        tensor_conn = [zeros_tensor(nconn, nv) for _ in props]
        props = (:pn, :sn, :bn)
        vec_cell = [zeros(Float64, nc) for _ in props]
        props = (:p_rec, :s_rec)
        vec_rec = [Dict{Float64, Vector{Float64}}() for _ in props]
        return new(tensor_cell..., tensor_conn..., vec_cell..., ρs, pvt, vec_rec...)
    end
end

struct Component
    a::Tensor   # Cells Component Accumulation
    r::Tensor   # Cells Component Residual
    function Component(nc::Int, nv::Int)
        props = (:a, :r)
        tensor_cell = [zeros_tensor(nc, nv) for _ in props]
        return new(tensor_cell...)
    end
end

function compute_ρ(phase::Phase)::Tensor
    phase.ρ .= phase.ρs ./ phase.b
end

function compute_λ(phase::Phase)::Tensor
    phase.λ .= phase.kr ./ (phase.μ .* phase.b)
end

function compute_b(phase::Phase)::Tensor where {T <: AbstractPVT}
    get_b(phase.pvt, phase.b, phase.p)
end

function compute_μ(phase::Phase)::Tensor where {T <: AbstractPVT}
    get_μ(phase.pvt, phase.μ, phase.p)
end


function compute_γ(phase::Phase, connlist::ConnList)::Tensor
    l_list, r_list = connlist.l, connlist.r
    ρ, γ = phase.ρ, phase.γ
    for i = 1:connlist.nconn
        l, r = l_list[i], r_list[i]
        γ[i] = β * g_ / gc * (ρ[l] + ρ[r]) / 2.
    end
    return γ
end

function compute_ΔΨ(phase::Phase, connlist::ConnList)::Tensor
    l_list, r_list, Δd = connlist.l, connlist.r, connlist.Δd
    ΔΨ, p, γ = phase.ΔΨ, phase.p, phase.γ
    for i = 1:connlist.nconn
        l, r = l_list[i], r_list[i]
        ΔΨ[i] = p[l] - p[r] - γ[i]*Δd[i]
    end
    return ΔΨ
end

function compute_f(phase::Phase, connlist::ConnList)::Tensor
    l_list, r_list, trans = connlist.l, connlist.r, connlist.trans
    ΔΨ, λ, f = phase.ΔΨ, phase.λ, phase.f
    for i = 1:connlist.nconn
        l, r = l_list[i], r_list[i]
        if ΔΨ[i] > 0.0
            f[i] = λ[l] * trans[i] * ΔΨ[i]
        else
            f[i] = λ[r] * trans[i] * ΔΨ[i]
        end
    end
    return f
end

function compute_a(
    phase::Phase,
    v::Vector{Float64},
    ϕ::Vector{Float64},
    dt::Float64,
)::Tensor
    return v .* ϕ .* (phase.s ./ phase.b - phase.sn ./ phase.bn) / dt / M
end


function set_phase_tn(phase::Phase, p::Vector{Float64}, s::Vector{Float64})
    phase.pn .= p
    phase.sn .= s
    get_b(phase.pvt, phase.bn, p)
end

function update_phase_tn(phase::Phase)::Phase
    for (x, y) in zip((:pn, :sn, :bn), (:p, :s, :b))
        getfield(phase, x) .= data(getfield(phase, y))
    end
    return phase
end

function update_phase(phase::Phase, connlist::ConnList)::Phase
    compute_b(phase)
    compute_μ(phase)
    compute_ρ(phase)
    compute_λ(phase)
    compute_γ(phase, connlist)
    compute_ΔΨ(phase, connlist)
    compute_f(phase, connlist)
    return phase
end


## Define Abstract State
abstract type AbstractFluid end

function update_fluid_tn(fluid::AbstractFluid)
    for phase in fluid.phases
        update_phase_tn(phase)
    end
end

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

function update_primary_variable(fluid::OWFluid, po::Vector{Float64}, sw::Vector{Float64})::OWFluid
    phases = fluid.phases
    phases.o.p .= param(po, 1, 2)
    phases.w.s .= param(sw, 2, 2)
    return fluid
end

function reset_primary_variable(fluid::OWFluid)
    phases = fluid.phases
    update_primary_variable(fluid, phases.o.pn, phases.w.sn)
end

function update_phases(fluid::OWFluid, connlist::ConnList)::OWFluid
    phases = fluid.phases
    phases.w.p .= phases.o.p
    phases.o.s .= 1 .- phases.w.s
    compute_kr(fluid)
    update_phase(fluid.phases.o, connlist)
    update_phase(fluid.phases.w, connlist)
    return fluid
end


end
