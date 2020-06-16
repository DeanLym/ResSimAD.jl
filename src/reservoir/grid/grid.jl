module Grid

using ..Global: α
using ..Rock: AbstractRock

export AbstractGrid,
    CartGrid, get_grid_index, set_cell_size, set_perm, set_poro, construct_conn

abstract type AbstractGrid end
# 1. Common fields of any concrete subtypes of AbstractGrid:
#       -- nc, inconn, outconn, v, ϕ, d, k
# 2. Common behaviors
#       -- input_{ϕ}

include("connlist.jl")

include("cartgrid.jl")

function set_cell_depth(grid::AbstractGrid, d::Vector{Float64})::AbstractGrid
    grid.d .= d
    return grid
end

function set_cell_depth(grid::AbstractGrid, d::Float64)::AbstractGrid
    set_cell_depth(grid, d*ones(grid.nc))
end

end
