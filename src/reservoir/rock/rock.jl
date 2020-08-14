module Rock

using Memento

const LOGGER = getlogger(@__MODULE__)
__init__() = Memento.register(LOGGER)

abstract type AbstractRock end

struct StandardRock <: AbstractRock
    nc::Int
    kx::Vector{Float64}
    ky::Vector{Float64}
    kz::Vector{Float64}
    ϕ::Vector{Float64}
    function StandardRock(nc::Int)
        props = (:kx, :ky, :kz, :ϕ)
        params = [ones(nc) for _ in props]
        return new(nc, params...)
    end
end

function set_perm(rock::StandardRock, kx::Float64, ky::Float64, kz::Float64)
    notice(LOGGER, "Setting constant permeability kx = $(round(kx, digits=3)) for all cells")
    notice(LOGGER, "Setting constant permeability ky = $(round(ky, digits=3)) for all cells")
    notice(LOGGER, "Setting constant permeability kz = $(round(kz, digits=3)) for all cells")
    I = ones(rock.nc)
    set_perm(rock, kx*I, ky*I, kz*I)
end

function set_perm(rock::StandardRock, kx::Vector{Float64}, ky::Vector{Float64}, kz::Vector{Float64})
    if any(kx .≤ 0.) || any(ky .≤ 0.) || any(kz .≤ 0.)
        error(LOGGER, "Negative value in permeability")
    end
    rock.kx .= kx
    rock.ky .= ky
    rock.kz .= kz
end

function set_poro(rock::StandardRock, ϕ::Float64)
    notice(LOGGER, "Setting constant porosity ϕ = $(round(ϕ, digits=3)) for all cells")
    I = ones(rock.nc)
    set_poro(rock, ϕ*I)
end

function set_poro(rock::StandardRock, ϕ::Vector{Float64})
    if any(ϕ .≤ 0.) error(LOGGER, "Negative value in porosity") end
    if any(ϕ .≥ 1.) error(LOGGER, "Porosity value larger than 1") end
    rock.ϕ .= ϕ
end

end
