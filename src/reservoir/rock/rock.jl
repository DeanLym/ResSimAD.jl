module Rock

abstract type AbstractRock end

struct StandardRock <: AbstractRock
    kx::Vector{Float64}
    ky::Vector{Float64}
    kz::Vector{Float64}
    ϕ::Vector{Float64}
    function StandardRock(nc::Int)
        props = (:kx, :ky, :kz, :ϕ)
        params = [ones(nc) for _ in props]
        return new(params...)
    end
end

function set_perm(rock::StandardRock, k::Vector{Float64})::StandardRock
    @assert all(k .> 0)
    rock.kx .= k
    rock.ky .= k
    rock.kz .= k
    return rock
end

function set_poro(grid::StandardRock, ϕ::Vector{Float64})::StandardRock
    @assert all(ϕ .> 0)
    grid.ϕ .= ϕ
    return grid
end

end
