module Grid

using ..Global: α

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

const Neighbor = Vector{Tuple{Vector{Int},Vector{Float64}}}

# Define struct for Conn List
struct ConnList
    numconn::Int
    l::Vector{Int}
    r::Vector{Int}
    trans::Vector{Float64}
    Δd::Vector{Float64} # Depth Diffrence
    function ConnList(numconn::Int)
        l = Vector{Int}(undef, numconn)
        r = Vector{Int}(undef, numconn)
        trans = Vector{Float64}(undef, numconn)
        Δd = Vector{Float64}(undef, numconn)
        return new(numconn, l, r, trans, Δd)
    end
end

function insert_conn(conn::ConnList, i::Int, l::Int, r::Int, trans::Float64, Δd::Float64)::ConnList
    conn.l[i] = l
    conn.r[i] = r
    conn.trans[i] = trans
    conn.Δd[i] = Δd
    return conn
end

# Define struct for Cartesian Grid
struct CartGrid <: AbstractStructGrid
    numcell::Int
    nx::Int
    ny::Int
    nz::Int
    connlist::ConnList

    dx::Vector{Float64}
    dy::Vector{Float64}
    dz::Vector{Float64}

    v::Vector{Float64} # Volume
    ϕ::Vector{Float64} # Porosity
    kx::Vector{Float64} # Permeability x
    ky::Vector{Float64} # Permeability y
    kz::Vector{Float64} # Permeability z
    d::Vector{Float64} # Depth
end

function CartGrid(nx::Int, ny::Int, nz::Int)::CartGrid
    numcell = nx*ny*nz
    numconn = (nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1)
    vecs = (:dx, :dy, :dz, :v, :ϕ, :kx, :ky, :kz, :d)
    params = [Vector{Float64}(undef, numcell) for v in vecs]
    # !format: off
    return CartGrid(numcell, nx, ny, nz, ConnList(numconn), params...)
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
    I = ones(Int, grid.numcell)
    set_cell_size(grid, dx * I, dy * I, dz * I)
end

function set_cell_depth(grid::AbstractGrid, d::Vector{Float64})::AbstractGrid
    grid.d .= d
    return grid
end

function set_cell_depth(grid::AbstractGrid, d::Float64)::AbstractGrid
    set_cell_depth(grid, d*ones(grid.numcell))
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

function set_perm(grid::AbstractGrid, k::Vector{Float64})::CartGrid
    @assert all(k .> 0)
    grid.kx .= k
    grid.ky .= k
    grid.kz .= k
    return grid
end

function set_perm(grid::AbstractGrid, k::Float64)::CartGrid
    @assert k > 0
    set_perm(grid, k*ones(grid.numcell))
end

function set_poro(grid::AbstractGrid, ϕ::Vector{Float64})::CartGrid
    @assert all(ϕ .> 0)
    grid.ϕ .= ϕ
    return grid
end

function set_poro(grid::AbstractGrid, ϕ::Float64)::CartGrid
    @assert ϕ > 0
    set_poro(grid, ϕ * ones(grid.numcell))
end

function construct_connlist(grid::CartGrid)::CartGrid
    count = 1
    dx, dy, dz = grid.dx, grid.dy, grid.dz
    kx, ky, kz = grid.kx, grid.ky, grid.kz #
    # X Direction
    for kk = 1:grid.nz
        for jj = 1:grid.ny
            for ii = 1:grid.nx-1
                l = get_grid_index(grid, ii, jj, kk)
                r = get_grid_index(grid, ii + 1, jj, kk)
                trans =
                    2 * α * dy[l] * dz[l] / ((dx[l] / kx[l]) + (dx[r] / kx[r]))
                Δd = grid.d[l] - grid.d[r]
                insert_conn(grid.connlist, count, l, r, trans, Δd)
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
                    2 * α * dx[l] * dz[l] / ((dy[l] / ky[l]) + (dy[r] / ky[r]))
                Δd = grid.d[l] - grid.d[r]
                insert_conn(grid.connlist, count, l, r, trans, Δd)
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
                    2 * α * dx[l] * dy[l] / ((dz[l] / kz[l]) + (dz[r] / kz[r]))
                Δd = grid.d[l] - grid.d[r]
                insert_conn(grid.connlist, count, l, r, trans, Δd)
                count += 1
            end
        end
    end
    return grid
end

end
