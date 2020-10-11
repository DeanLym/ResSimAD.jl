module Rock

using Memento

const LOGGER = getlogger(@__MODULE__)
__init__() = Memento.register(LOGGER)

abstract type AbstractRock end

struct StandardRock <: AbstractRock
    nc::Int
    kx::Vector{Float64} # x permeability 
    ky::Vector{Float64} # y permeability 
    kz::Vector{Float64} # z permeability 
    ϕ::Vector{Float64} # porosity
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

function set_poro(rock::T, ϕ::Float64) where{T <: AbstractRock}
    notice(LOGGER, "Setting constant porosity ϕ = $(round(ϕ, digits=3)) for all cells")
    I = ones(rock.nc)
    set_poro(rock, ϕ*I)
end

function set_poro(rock::T, ϕ::Vector{Float64}) where{T <: AbstractRock}
    if any(ϕ .≤ 0.) error(LOGGER, "Negative value in porosity") end
    if any(ϕ .≥ 1.) error(LOGGER, "Porosity value larger than 1") end
    rock.ϕ .= ϕ
end

struct TransRock <: AbstractRock
    nc::Int
    tranx::Vector{Float64} # x Transmissibility
    trany::Vector{Float64} # y Transmissibility
    tranz::Vector{Float64} # z Transmissibility
    ϕ::Vector{Float64} # porosity
    function TransRock(nc::Int)
        props = (:tranx, :trany, :tranz, :ϕ)
        params = [ones(nc) for _ in props]
        return new(nc, params...)
    end 
end

function set_trans(rock::TransRock, tranx::Float64, trany::Float64, tranz::Float64)
    notice(LOGGER, "Setting constant transmissibility tranx = $(round(tranx, digits=3)) for all cells")
    notice(LOGGER, "Setting constant transmissibility trany = $(round(trany, digits=3)) for all cells")
    notice(LOGGER, "Setting constant transmissibility tranz = $(round(tranz, digits=3)) for all cells")
    I = ones(rock.nc)
    set_perm(rock, tranx*I, trany*I, tranz*I)
end

function set_trans(rock::TransRock, tranx::Vector{Float64}, trany::Vector{Float64}, tranz::Vector{Float64})
    if any(tranx .< 0.) || any(trany .< 0.) || any(tranz .< 0.)
        error(LOGGER, "Negative value in transmissibility")
    end
    rock.tranx .= tranx
    rock.trany .= trany
    rock.tranz .= tranz
end

end
