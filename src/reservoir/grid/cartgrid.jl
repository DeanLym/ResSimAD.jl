abstract type AbstractStructGrid <: AbstractGrid end
# 1. Common fields of any concreate subtypes of AbstractStructGrid
#       -- nx, ny, nz, k
# 2. Common behaviors
#       -- input_{k}
#       -- getgridindex(i,j,k -> ind)

# Define struct for Cartesian Grid
struct CartGrid <: AbstractStructGrid
    nc::Int
    nx::Int
    ny::Int
    nz::Int

    dx::Vector{Float64}
    dy::Vector{Float64}
    dz::Vector{Float64}

    v::Vector{Float64} # Volume
    d::Vector{Float64} # Depth

    connlist::ConnList
    neighbors::Vector{Vector{Int64}}
end

function CartGrid(nx::Int, ny::Int, nz::Int)::CartGrid
    nc = nx*ny*nz
    vecs = (:dx, :dy, :dz, :v, :d)
    params = [Vector{Float64}(undef, nc) for v in vecs]
    neighbors = Vector{Vector{Int}}(undef, nc)
    # !format: off
    return CartGrid(nc, nx, ny, nz, params..., ConnList(), neighbors)
    # !format: on
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
)::CartGrid
    @assert dx > 0 && dy > 0 && dz > 0
    I = ones(Int, grid.nc)
    set_cell_size(grid, dx * I, dy * I, dz * I)
end


function set_cell_size(
    grid::CartGrid,
    dx::Vector{Float64},
    dy::Vector{Float64},
    dz::Vector{Float64},
)::CartGrid
    @assert all(dx .> 0) && all(dy .> 0) && all(dz .> 0)
    grid.dx .= dx
    grid.dy .= dy
    grid.dz .= dz
    # Compute volume
    grid.v .= grid.dx .* grid.dy .* grid.dz
    return grid
end


function construct_connlist(grid::CartGrid, rock::AbstractRock)::ConnList
    connlist = grid.connlist
    dx, dy, dz = grid.dx, grid.dy, grid.dz
    kx, ky, kz = rock.kx, rock.ky, rock.kz #
    # X Direction
    for kk = 1:grid.nz
        for jj = 1:grid.ny
            for ii = 1:grid.nx-1
                l = get_grid_index(grid, ii, jj, kk)
                r = get_grid_index(grid, ii + 1, jj, kk)
                trans =
                    2 * α * dy[l] * dz[l] / ((dx[l] / kx[l]) + (dx[r] / kx[r]))
                Δd = grid.d[l] - grid.d[r]
                add_conn(connlist, l, r, trans, Δd)
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
                    2 * α * dx[l] * dz[l] / ((dy[l] / ky[l]) + (dy[r] / ky[r]))
                Δd = grid.d[l] - grid.d[r]
                add_conn(connlist, l, r, trans, Δd)
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
                    2 * α * dx[l] * dy[l] / ((dz[l] / kz[l]) + (dz[r] / kz[r]))
                Δd = grid.d[l] - grid.d[r]
                add_conn(connlist, l, r, trans, Δd)
            end
        end
    end
    connlist.nconn = length(connlist.trans)
    construct_neighbors(grid)
    return connlist
end
