module AutoDiff

using StaticArrays
import Base
import Calculus
using SpecialFunctions

export param, zeros_tensor, Tensor

const SVec = SVector{Nv, Float64} where {Nv}
const SMat = SMatrix{Nv, Nc, Float64, L} where {Nv, Nc, L}

struct SGrad{Nv}
    ic::Int
    ∇::SVec{Nv}
end

struct DGrad{Nv, L}
    ic::SVector{2, Int}
    ∇::SMat{Nv, 2, L}
end

struct MGrad{Nv, Nc, L}
    ic::Int
    ∇::SVec{Nv}
    icn::SVector{Nc, Int}
    ∇n::SMat{Nv, Nc, L}
end

MGrad(ic::Int, nv::Int) =
    MGrad{nv,0,0}(ic, zeros(SVec{nv}), SVector{0,Int}(), SMat{nv,0,0}())

const Grad{Nv, Nc, L} = Union{SGrad{Nv}, DGrad{Nv, L}, MGrad{Nv, Nc}} where {Nv, Nc, L}

mutable struct Variable{T<:Grad} <: Number
    val::Float64
    grad::T
    Variable{T}(val::Float64, grad::T) where {T <: Grad} = new(val, grad)
end

const Tensor = Union{Vector{Float64}, Vector{Variable}, Vector{Variable{T}}}  where{T<:Grad}

Variable(val::Float64, grad::SGrad{Nv}) where {Nv} =
    Variable{SGrad{Nv}}(val, grad)

Variable(val::Float64, grad::DGrad{Nv,L}) where {Nv,L} =
    Variable{DGrad{Nv,L}}(val, grad)

Variable(val::Float64, grad::MGrad{Nv,Nc, L}) where {Nv,Nc,L} =
    Variable{MGrad{Nv,Nc,L}}(val, grad)

function Variable(val::Float64, ic::Int, grad::SVector{Nv, Float64}) where {Nv}
    Variable(val, SGrad{Nv}(ic, grad))
end


function param(x::Vector{Float64}, iv::Int, nv::Int)
    grad = SVector{nv, Float64}([if i==iv 1.0 else 0.0 end for i=1:nv])
    n = length(x)
    return [Variable(x[i], i, grad) for i=1:n]
end

function zeros_tensor(nc::Int, nv::Int)
    grad = zeros(SVector{nv, Float64})
    vec = Vector{Variable}([Variable(0.0, i, grad) for i=1:nc])
    return vec
end

function zeros_tensor(ind::Vector{Int}, nv::Int)
    grad = zeros(SVector{nv, Float64})
    vec = Vector{Variable}([Variable(0.0, i, grad) for i in ind])
    return vec
end

function ones_tensor(nc::Int, nv::Int)
    grad = zeros(SVector{nv, Float64})
    vec = Vector{Variable}([Variable(1.0, i, grad) for i=1:nc])
    return vec
end

function ones_tensor(ind::Vector{Int}, nv::Int)
    grad = zeros(SVector{nv, Float64})
    vec = Vector{Variable}([Variable(1.0, i, grad) for i in ind])
    return vec
end

function grad(x::Variable{MGrad{Nv, Nc, L}}) where {Nv, Nc, L}
    g = x.grad
    return vcat(SA[g.ic], g.icn), hcat(g.∇, g.∇n)
end

function grad(x::Variable{DGrad{Nv, L}}) where {Nv, L}
    return x.grad.ic, x.grad.∇
end

function grad(x::Variable{SGrad{Nv}}) where {Nv}
    return SA[x.grad.ic], SMat{Nv, 1, Nv}(x.grad.∇)
end

function data(t::Tensor)::Vector{Float64}
    return [x.val for x in t]
end

## Operator overload
# SGrad and Real
function Base.:-(z::SGrad{Nv})::SGrad{Nv} where {Nv}
    return SGrad{Nv}(z.ic, -z.∇)
end

function Base.:*(z::SGrad{Nv}, w::Real)::SGrad{Nv} where {Nv}
    return SGrad{Nv}(z.ic, z.∇*w)
end

function Base.:*(w::Real, z::SGrad{Nv})::SGrad{Nv} where {Nv}
    return SGrad{Nv}(z.ic, w*z.∇)
end

function Base.:/(z::SGrad{Nv}, w::Real)::SGrad{Nv} where {Nv}
    return SGrad{Nv}(z.ic, z.∇/w)
end

function Base.:/(w::Real, z::SGrad{Nv})::SGrad{Nv} where {Nv}
    return SGrad{Nv}(z.ic, w/z.∇)
end

# SGrad and SGrad
function Base.:+(z::SGrad{Nv}, w::SGrad{Nv}) where {Nv}
    if z.ic == w.ic
        return SGrad{Nv}(z.ic, z.∇ + w.∇)
    elseif z.ic < w.ic
        return DGrad{Nv, 2*Nv}(SA[z.ic, w.ic], hcat(z.∇, w.∇))
    else
        return DGrad{Nv, 2*Nv}(SA[w.ic, z.ic], hcat(w.∇, z.∇))
    end
end

function Base.:-(z::SGrad{Nv}, w::SGrad{Nv}) where {Nv}
    if z.ic == w.ic
        return SGrad{Nv}(z.ic, z.∇ - w.∇)
    elseif z.ic < w.ic
        return DGrad{Nv, 2*Nv}(SA[z.ic, w.ic], hcat(z.∇, -w.∇))
    else
        return DGrad{Nv, 2*Nv}(SA[w.ic, z.ic], hcat(-w.∇, z.∇))
    end
end

# DGrad Real
function Base.:-(z::DGrad{Nv, L})::DGrad{Nv, L} where {Nv, L}
    return DGrad{Nv, L}(z.ic, -z.∇)
end

function Base.:*(z::DGrad{Nv, L}, w::Real)::DGrad{Nv, L} where {Nv, L}
    return DGrad{Nv, L}(z.ic, z.∇*w)
end

function Base.:*(w::Real, z::DGrad{Nv, L})::DGrad{Nv, L} where {Nv, L}
    return DGrad{Nv, L}(z.ic, w*z.∇)
end

function Base.:/(z::DGrad{Nv, L}, w::Real)::DGrad{Nv, L} where {Nv, L}
    return DGrad{Nv, L}(z.ic, z.∇/w)
end

function Base.:/(w::Real, z::DGrad{Nv, L})::DGrad{Nv, L} where {Nv, L}
    return DGrad{Nv, L}(z.ic, w/z.∇)
end

# SGrad and DGrad
function Base.:+(z::SGrad{Nv}, w::DGrad{Nv, L})::DGrad{Nv, L} where {Nv, L}
    g = zeros(SVector{Nv, Float64})
    if z.ic == w.ic[1]
        return DGrad{Nv, L}(w.ic, hcat(z.∇, g) + w.∇)
    else
        return DGrad{Nv, L}(w.ic, hcat(g, z.∇) + w.∇)
    end
end

function Base.:+(w::DGrad{Nv, L}, z::SGrad{Nv})::DGrad{Nv, L} where {Nv, L}
    g = zeros(SVector{Nv, Float64})
    if z.ic == w.ic[1]
        return DGrad{Nv, L}(w.ic, hcat(z.∇, g) + w.∇)
    else
        return DGrad{Nv, L}(w.ic, hcat(g, z.∇) + w.∇)
    end
end

function Base.:-(z::SGrad{Nv}, w::DGrad{Nv, L})::DGrad{Nv, L} where {Nv, L}
    g = zeros(SVector{Nv, Float64})
    if z.ic == w.ic[1]
        return DGrad{Nv, L}(w.ic, hcat(z.∇, g) - w.∇)
    else
        return DGrad{Nv, L}(w.ic, hcat(g, z.∇) - w.∇)
    end
end

function Base.:-(w::DGrad{Nv, L}, z::SGrad{Nv})::DGrad{Nv, L} where {Nv, L}
    g = zeros(SVector{Nv, Float64})
    if z.ic == w.ic[1]
        return DGrad{Nv, L}(w.ic, hcat(-z.∇, g) + w.∇)
    else
        return DGrad{Nv, L}(w.ic, hcat(g, -z.∇) + w.∇)
    end
end

##
# DGrad and DGrad
function Base.:+(z::DGrad{Nv,L}, w::DGrad{Nv,L}) where {Nv,L}
    zic, wic = z.ic, w.ic
    z∇, w∇ = z.∇, w.∇
    if zic == wic #zic and wic are sorted
        return DGrad{Nv, L}(zic, z∇ + w∇)
    elseif zic[2] == wic[1]
        return MGrad{Nv,2,2 * Nv}(
            zic[2],
            z∇[:, 2] + w∇[:, 1],
            SA[zic[1], wic[2]],
            hcat(z∇[:, 1], w∇[:, 2]),
        )
    elseif zic[1] == wic[2]
        return MGrad{Nv,2,2 * Nv}(
            zic[1],
            z∇[:, 1] + w∇[:, 2],
            SA[zic[2], wic[1]],
            hcat(z∇[:, 2], w∇[:, 1]),
        )
    elseif zic[1] == wic[1]
        return MGrad{Nv,2,2 * Nv}(
            zic[1],
            z∇[:, 1] + w∇[:, 1],
            SA[zic[2], wic[2]],
            hcat(z∇[:, 2], w∇[:, 2]),
        )
    else
        return MGrad{Nv,2,2 * Nv}(
            zic[2],
            z∇[:, 2] + w∇[:, 2],
            SA[zic[1], wic[1]],
            hcat(z∇[:, 1], w∇[:, 1]),
        )
    end
end

function Base.:-(z::DGrad{Nv,L}, w::DGrad{Nv,L}) where {Nv,L}
    zic, wic = z.ic, w.ic
    z∇, w∇ = z.∇, w.∇
    if zic == wic #zic and wic are sorted
        return DGrad{Nv, L}(zic, z∇ - w∇)
    elseif zic[2] == wic[1]
        return MGrad{Nv,2,2 * Nv}(
            zic[2],
            z∇[:, 2] - w∇[:, 1],
            SA[zic[1], wic[2]],
            hcat(z∇[:, 1], -w∇[:, 2]),
        )
    elseif zic[1] == wic[2]
        return MGrad{Nv,2,2 * Nv}(
            zic[1],
            z∇[:, 1] - w∇[:, 2],
            SA[zic[2], wic[1]],
            hcat(z∇[:, 2], -w∇[:, 1]),
        )
    elseif zic[1] == wic[1]
        return MGrad{Nv,2,2 * Nv}(
            zic[1],
            z∇[:, 1] - w∇[:, 1],
            SA[zic[2], wic[2]],
            hcat(z∇[:, 2], -w∇[:, 2]),
        )
    else
        return MGrad{Nv,2,2 * Nv}(
            zic[2],
            z∇[:, 2] - w∇[:, 2],
            SA[zic[1], wic[1]],
            hcat(z∇[:, 1], -w∇[:, 1]),
        )
    end
end

## MGrad and SGrad
function Base.:+(z::MGrad{Nv, Nc, Lz}, w::SGrad{Nv}) where {Nv, Nc, Lz}
    return MGrad{Nv,Nc, Lz}(z.ic,z.∇ + w.∇, z.icn, z.∇n)
end

function Base.:+(w::SGrad{Nv}, z::MGrad{Nv, Nc, Lz}) where {Nv, Nc, Lz}
    return MGrad{Nv,Nc, Lz}(z.ic,z.∇ + w.∇, z.icn, z.∇n)
end

function Base.:-(z::MGrad{Nv, Nc, Lz}, w::SGrad{Nv}) where {Nv, Nc, Lz}
    return MGrad{Nv,Nc, Lz}(z.ic,z.∇ - w.∇, z.icn, z.∇n)
end

function Base.:-(w::SGrad{Nv}, z::MGrad{Nv, Nc, Lz}) where {Nv, Nc, Lz}
    return MGrad{Nv,Nc, Lz}(z.ic,w.∇ - z.∇, z.icn, z.∇n)
end

##  MGrad and DGrad
function Base.:+(z::MGrad{Nv, Nc, Lz}, w::DGrad{Nv,Lw}) where {Nv, Nc, Lz,Lw}
    zic, wic = z.ic, w.ic
    w∇ = w.∇
    if zic == wic[1]
        return MGrad{Nv,Nc+1, Nv*(Nc+1)}(
            zic,
            z.∇ + w∇[:, 1],
            vcat(z.icn, SA[wic[2]]),
            hcat(z.∇n, w∇[:, 2]),
        )
    else
        return MGrad{Nv,Nc+1, Nv*(Nc+1)}(
            zic,
            z.∇ + w∇[:, 2],
            vcat(z.icn, SA[wic[1]]),
            hcat(z.∇n, w∇[:, 1]),
        )
    end
end

function Base.:+(w::DGrad{Nv,Lw}, z::MGrad{Nv, Nc, Lz}) where {Nv, Nc, Lz,Lw}
    zic, wic = z.ic, w.ic
    w∇ = w.∇
    if zic == wic[1]
        return MGrad{Nv,Nc+1, Nv*(Nc+1)}(
            zic,
            z.∇ + w∇[:, 1],
            vcat(z.icn, SA[wic[2]]),
            hcat(z.∇n, w∇[:, 2]),
        )
    else
        return MGrad{Nv,Nc+1, Nv*(Nc+1)}(
            zic,
            z.∇ + w∇[:, 2],
            vcat(z.icn, SA[wic[1]]),
            hcat(z.∇n, w∇[:, 1]),
        )
    end
end

function Base.:-(z::MGrad{Nv, Nc, Lz}, w::DGrad{Nv,Lw}) where {Nv, Nc, Lz,Lw}
    zic, wic = z.ic, w.ic
    w∇ = w.∇
    if zic == wic[1]
        return MGrad{Nv,Nc+1, Nv*(Nc+1)}(
            zic,
            z.∇ - w∇[:, 1],
            vcat(z.icn, SA[wic[2]]),
            hcat(z.∇n, -w∇[:, 2]),
        )
    else
        return MGrad{Nv,Nc+1, Nv*(Nc+1)}(
            zic,
            z.∇ - w∇[:, 2],
            vcat(z.icn, SA[wic[1]]),
            hcat(z.∇n, -w∇[:, 1]),
        )
    end
end

function Base.:-(w::DGrad{Nv,Lw}, z::MGrad{Nv, Nc, Lz}) where {Nv, Nc, Lz,Lw}
    zic, wic = z.ic, w.ic
    w∇ = w.∇
    if zic == wic[1]
        return MGrad{Nv,Nc+1, Nv*(Nc+1)}(
            zic,
            w∇[:, 1] - z.∇,
            vcat(z.icn, SA[wic[2]]),
            hcat(-z.∇n, w∇[:, 2]),
        )
    else
        return MGrad{Nv,Nc+1, Nv*(Nc+1)}(
            zic,
            w∇[:, 2] - z.∇,
            vcat(z.icn, SA[wic[1]]),
            hcat(-z.∇n, w∇[:, 1]),
        )
    end
end

## Auto Diff Rule

Base.:(==)(z::Variable, w::Variable)   = z.val == w.val
Base.:(==)(z::Variable, x::Real)   = z.val == x
Base.:(==)(x::Real, z::Variable)   = z.val == x

Base.isequal(z::Variable, w::Variable)  = isequal(z.val,w.val) && isequal(z.grad, w.grad)
Base.isequal(z::Variable, x::Real)   = isequal(z.val, x) && isempty(z.grad)
Base.isequal(x::Real, z::Variable)   = isequal(z, x)

Base.isless(z::Variable,w::Variable)   = z.val < w.val
Base.isless(z::Real,w::Variable)   = z < w.val
Base.isless(z::Variable,w::Real)   = z.val < w

Base.floor(z::Variable)  = floor(z.val)
Base.ceil(z::Variable)   = ceil(z.val)
Base.trunc(z::Variable)  = trunc(z.val)
Base.round(z::Variable)  = round(z.val)
Base.floor(::Type{T}, z::Variable) where {N, T<:Real} = floor(T, z.val)
Base.trunc(::Type{T}, z::Variable) where {N, T<:Real} = trunc(T, z.val)
Base.ceil( ::Type{T}, z::Variable) where {N, T<:Real} = ceil( T, z.val)
Base.round(::Type{T}, z::Variable) where {N, T<:Real} = round(T, z.val)

Base.abs(z::Variable) = z ≥ 0 ? z : -z

Base.:+(z::Variable, w::Variable)   = Variable(z.val + w.val, z.grad +w.grad)
Base.:+(z::Real, w::Variable)   = Variable(z+w.val, w.grad)
Base.:+(z::Variable, w::Real)   = Variable(z.val+w, z.grad)

Base.:-(z::Variable)   = Variable(-z.val, -z.grad)
Base.:-(z::Variable, w::Variable)   = Variable(z.val - w.val, z.grad - w.grad)
Base.:-(z::Real, w::Variable)   = Variable(z-w.val, -w.grad)
Base.:-(z::Variable, w::Real)   = Variable(z.val-w, z.grad)

Base.:*(z::Variable, w::Variable)   = Variable(z.val * w.val, z.grad*w.val+z.val*w.grad)
Base.:*(x::Real, z::Variable)   = Variable(x*z.val, x*z.grad)
Base.:*(z::Variable, x::Number)   = Variable(x*z.val, x*z.grad)

Base.:/(z::Variable, w::Variable)   = Variable(z.val/w.val, (z.grad*w.val-z.val*w.grad)/(w.val*w.val))
Base.:/(z::Real, w::Variable)   = Variable(z/w.val, -z*w.grad/w.val^2)
Base.:/(z::Variable, x::Real)   = Variable(z.val/x, z.grad/x)


Base.:^(z::Variable, n::Integer)  = Variable(z.val^n, z.grad*n*z.val^(n-1))
Base.:^(z::Variable, n::Rational)   = Variable(z.val^n, z.grad*n*z.val^(n-1))

Base.:^(z::Variable, n::Real)   = Variable(z.val^n, z.grad*n*z.val^(n-1))

function Base.:^(z::Variable, w::Variable)
    val = Base.:^(z.val, w.val)
    grad = z.grad * w.val * Base.:^(z.val, w.val - 1) +
         w.grad * Base.:^(z.val, w.val) * log(z.val)
    return Variable(val, grad)
end

Base.mod(z::Variable, n::Real)   = Variable(mod(z.val, n), z.grad)

Base.inv(z::Variable)  = Variable(inv(z.val),-z.grad/z.val^2)


for (funsym, exp) in Calculus.symbolic_derivatives_1arg()
    funsym == :exp && continue
    funsym == :abs2 && continue
    funsym == :inv && continue
    if isdefined(SpecialFunctions, funsym)
        @eval function SpecialFunctions.$(funsym)(z::Variable)
            x = z.val
            xp = z.grad
            Variable($(funsym)(x),xp*$exp)
        end
    elseif isdefined(Base, funsym)
        @eval function Base.$(funsym)(z::Variable)
            x = z.val
            xp = z.grad
            Variable($(funsym)(x),xp*$exp)
        end
    end
end

# only need to compute exp/cis once
Base.exp(z::Variable)  = (expval = exp(z.val); Variable(expval, z.grad*expval))
Base.exp10(x::Variable)  = (y = exp10(x.val); Variable(y, y * log(10) * x.grad))

Base.sinpi(z::Variable)  = Variable(sinpi(z.val), z.grad*cospi(z.val)*π)
Base.cospi(z::Variable)  = Variable(cospi(z.val), -z.grad*sinpi(z.val)*π)

end
