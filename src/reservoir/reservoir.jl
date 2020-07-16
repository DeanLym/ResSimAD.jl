module Reservoir

using ..AutoDiff: value
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
        phase.p_rec[round(t, digits=6)] = value(phase.p)
        phase.s_rec[round(t, digits=6)] = value(phase.s)
    end
end



end
