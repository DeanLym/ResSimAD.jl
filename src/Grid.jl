module Grid

export AbstractGrid,
    CartGrid, get_grid_index, set_cell_size, set_perm, set_poro, construct_conn

abstract type AbstractGrid end
# 1. Common fields of any concrete subtypes of AbstractGrid:
#       -- numcell, inconn, outconn, v, ϕ, d, k
# 2. Common behaviors
#       -- input_{ϕ}
abstract type AbstractStructGrid <: AbstractGrid end
# 1. Common fields of any concreate subtypes of AbstractStructGrid
#       -- nx, ny, nz, k
# 2. Common behaviors
#       -- input_{k}
#       -- getgridindex(i,j,k -> ind)
const α = 0.001127
const Neighbor = Vector{Tuple{Vector{Int},Vector{Float64}}}

# Define struct for Conn List
struct ConnList
    numconn::Int
    l::Vector{Int}
    r::Vector{Int}
    trans::Vector{Float64}
    function ConnList(numconn::Int)
        l = Vector{Int}(undef, numconn)
        r = Vector{Int}(undef, numconn)
        trans = Vector{Float64}(undef, numconn)
        return new(numconn, l, r, trans)
    end
end

function insert_conn(conn::ConnList, ind::Int, l::Int, r::Int, trans::Float64)
    conn.l[ind] = l
    conn.r[ind] = r
    conn.trans[ind] = trans
end


# Define struct for Cartesian Grid
struct CartGrid <: AbstractStructGrid
    numcell::Int
    nx::Int
    ny::Int
    nz::Int

    dx::Vector{Float64}
    dy::Vector{Float64}
    dz::Vector{Float64}

    v::Vector{Float64} # Volume
    ϕ::Vector{Float64} # Porosity
    k::Vector{Float64} # Permeability
    d::Vector{Float64} # Depth

    connlist::ConnList
    neighbors::Neighbor # Neighbors for each cell
    # inconn::Array{Array{Tuple{Int, Float64}, 1}, 1}  # Inbound connection for each cell
    # outconn::Array{Array{Tuple{Int, Float64}, 1}, 1} # Outbound connection for each cell

    function CartGrid(nx::Int, ny::Int, nz::Int)::CartGrid
        numcell = nx * ny * nz

        dx = Vector{Float64}(undef, numcell)
        dy = Vector{Float64}(undef, numcell)
        dz = Vector{Float64}(undef, numcell)

        v = Vector{Float64}(undef, numcell)
        ϕ = Vector{Float64}(undef, numcell)
        k = Vector{Float64}(undef, numcell)
        d = Vector{Float64}(undef, numcell)

        numconn = (nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1)
        connlist = ConnList(numconn)

        neighbors = Neighbor(undef, numcell)

        return new(
            numcell,
            nx,
            ny,
            nz,
            dx,
            dy,
            dz,
            v,
            ϕ,
            k,
            d,
            connlist,
            neighbors,
        )
    end
end

function get_grid_index(grid::AbstractStructGrid, i::Int, j::Int, k::Int)::Int
    @assert 0 < i <= grid.nx && 0 < j <= grid.ny && 0 < k <= grid.nz
    return grid.nx * grid.ny * (k - 1) + grid.nx * (j - 1) + i
end

function set_cell_size(
    grid::CartGrid,
    dx::Float64,
    dy::Float64,
    dz::Float64,
)::Nothing
    @assert dx > 0 && dy > 0 && dz > 0
    grid.dx .= dx
    grid.dy .= dy
    grid.dz .= dz
    # Compute volume
    grid.v .= grid.dx .* grid.dy .* grid.dz
    return nothing
end

function set_cell_size(
    grid::CartGrid,
    dx::Vector{Float64},
    dy::Vector{Float64},
    dz::Vector{Float64},
)::Nothing
    @assert all(dx .> 0) && all(dy .> 0) && all(dz .> 0)
    grid.dx .= dx
    grid.dy .= dy
    grid.dz .= dz
    # Compute volume
    grid.v .= grid.dx .* grid.dy .* grid.dz
    return nothing
end

function set_perm(grid::AbstractGrid, k::Vector{Float64})::Nothing
    @assert all(k .> 0)
    grid.k .= k
    return nothing
end

function set_perm(grid::AbstractGrid, k::Float64)::Nothing
    @assert k > 0
    grid.k .= k
    return nothing
end

function set_poro(grid::AbstractGrid, ϕ::Vector{Float64})::Nothing
    @assert all(ϕ .> 0)
    grid.ϕ .= ϕ
    return nothing
end

function set_poro(grid::AbstractGrid, ϕ::Float64)::Nothing
    @assert ϕ > 0
    grid.ϕ .= ϕ
    return nothing
end

function construct_conn(grid::CartGrid)::Nothing
    count = 1
    dx, dy, dz = grid.dx, grid.dy, grid.dz
    k = grid.k #
    # X Direction
    for kk = 1:grid.nz
        for jj = 1:grid.ny
            for ii = 1:grid.nx-1
                l = get_grid_index(grid, ii, jj, kk)
                r = get_grid_index(grid, ii + 1, jj, kk)
                trans =
                    2 * α * dy[l] * dz[l] / ((dx[l] / k[l]) + (dx[r] / k[r]))
                insert_conn(grid.connlist, count, l, r, trans)
                count += 1
            end
        end
    end
    # Y Direction
    for kk = 1:grid.nz
        for ii = 1:grid.nx
            for jj = 1:grid.ny-1
                l = get_grid_index(grid, ii, jj, kk)
                r = get_grid_index(grid, ii, jj + 1, kk)
                trans =
                    2 * α * dx[l] * dz[l] / ((dy[l] / k[l]) + (dy[r] / k[r]))
                insert_conn(grid.connlist, count, l, r, trans)
                count += 1
            end
        end
    end
    # Z Direction
    for jj = 1:grid.ny
        for ii = 1:grid.nx
            for kk = 1:grid.nz-1
                l = get_grid_index(grid, ii, jj, kk)
                r = get_grid_index(grid, ii, jj, kk + 1)
                trans =
                    2 * α * dx[l] * dy[l] / ((dz[l] / k[l]) + (dz[r] / k[r]))
                insert_conn(grid.connlist, count, l, r, trans)
                count += 1
            end
        end
    end
    # Initialize neighbors
    for ind = 1:grid.numcell
        grid.neighbors[ind] = (Int[], Float64[])
    end
    # Find neighbors for each cell
    for ind = 1:grid.connlist.numconn
        push!(grid.neighbors[grid.connlist.l[ind]][1], grid.connlist.r[ind])
        push!(grid.neighbors[grid.connlist.l[ind]][2], grid.connlist.trans[ind])
        push!(grid.neighbors[grid.connlist.r[ind]][1], grid.connlist.l[ind])
        push!(grid.neighbors[grid.connlist.r[ind]][2], grid.connlist.trans[ind])
    end
    return nothing
end

end
