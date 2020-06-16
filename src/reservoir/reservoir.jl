module Reservoir

using ..AutoDiff: data
using ..Rock: AbstractRock, StandardRock
using ..Grid: AbstractGrid, CartGrid, ConnList
using ..Fluid: AbstractFluid, Phase

abstract type AbstractReservoir end

struct StandardReservoir <: AbstractReservoir
    grid::AbstractGrid
    rock::AbstractRock
    fluid::AbstractFluid
end

function save_fluid_results(reservoir::AbstractReservoir, t::Float64)
    for phase in reservoir.fluid.phases
        push!(phase.p_rec, data(phase.p))
        push!(phase.s_rec, data(phase.s))
    end
end

end
