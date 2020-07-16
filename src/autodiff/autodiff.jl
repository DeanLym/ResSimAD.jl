module AutoDiff

using StaticArrays
import Base
const SVec = SVector{Nv, Float64} where {Nv}
const SMat = SMatrix{Nv, Nn, Float64, L} where {Nv, Nn, L}

struct ADScaler{Nv, Nn, L} <: Number
    v::Float64
    i::SVector{Nn, Int}
    ∇::SMat{Nv, Nn, L}
end

const ADVector = Vector{T} where T <: ADScaler

## Constructors
function advector(v::Vector{Float64}, g::SMat{Nv, 1, L}) where{Nv, L}
    return [ADScaler{Nv, 1, L}(v[i], SA[i], g) for i in eachindex(v)]
end

function advector(v::Vector{Float64}, iv::Int, nv::Int)
    g = SMat{nv, 1, nv}([if i==iv 1.0 else 0.0 end for i=1:nv])
    return advector(v, g)
end

function advector(n::Int, nn::Int, nv::Int)
    # n - vector length
    # nn - number of neighbors
    # nv - number of variables
    return Vector{ADScaler{nv, nn, nv*nn}}(undef, n)
end

## Utility function
function value(x::Tx) where {Tx <: ADScaler}
    return x.v
end

function value(x::Tv) where {Tv <: ADVector}
    return @. value(x)
end

function grad(x::Tx, in::Int, iv::Int) where {Tx <: ADScaler}
    # Get gradient of x w.r.p to the [iv]-th variable of the [in]-th neighbor
    return x.∇[iv, in]
end

function grad(v::Tv, in::Int, iv::Int) where {Tv <: ADVector}
    return @. grad(v, in, iv)
end

function index(x::Tx, in::Int) where {Tx <: ADScaler}
    # Get cell index for the [in]-th neighbor
    return x.i[in]
end

function index(v::Tv, in::Int) where {Tv <: ADVector}
    return @. index(v, in)
end

function adscaler(v::Float64, i::Int64, g::SMat{Nv, 1, Nv}) where {Nv}
    return ADScaler{Nv, 1, Nv}(v, SA[i], g)
end

function set_primary_variable(x::Vector{ADScaler{Nv, 1, Nv}}, v::Vector{Float64}, iv::Int) where {Nv}
    nc = length(v)
    g = [SMat{Nv, 1, Nv}([if i==iv 1.0 else 0.0 end for i=1:Nv])]
    @. x = adscaler(v, 1:nc, g)
end


## Operator Overload
function Base.:-(z::Tv) where {Tv <: ADScaler}
    return Tv(-z.v, z.i, -z.∇)
end

function Base.:+(z::Tv, w::T)::Tv where {Tv<:ADScaler, T<:Real}
    return Tv(w+z.v, z.i, z.∇)
end

function Base.:+(w::T, z::Tv)::Tv where {Tv<:ADScaler, T<:Real}
    return Tv(w+z.v, z.i, z.∇)
end

function Base.:-(z::Tv, w::T)::Tv where {Tv<:ADScaler, T<:Real}
    return Tv(z.v - w, z.i, z.∇)
end

function Base.:-(w::T, z::Tv)::Tv where {Tv<:ADScaler, T<:Real}
    return Tv(w - z.v, z.i, -z.∇)
end

function Base.:*(z::Tv, w::T)::Tv where {Tv<:ADScaler, T<:Real}
    return Tv(w*z.v, z.i, w*z.∇)
end

function Base.:*(w::T, z::Tv)::Tv where {Tv<:ADScaler, T<:Real}
    return Tv(w*z.v, z.i, w*z.∇)
end

function Base.:/(z::Tv, w::T)::Tv where {Tv<:ADScaler, T<:Real}
    return Tv(z.v/w, z.i, z.∇/w)
end

function Base.:/(w::T, z::Tv)::Tv where {Tv<:ADScaler, T<:Real}
    return Tv(w/z.v, z.i, -w*z.∇/z.v^2)
end

function Base.:^(z::Tv, n::Integer) where {Tv<:ADScaler}
    return Tv(z.v^n, z.i, z.∇*n*z.v^(n-1))
end

function Base.:^(z::Tv, n::Rational) where {Tv<:ADScaler}
    return Tv(z.v^n, z.i, z.∇*n*z.v^(n-1))
end

function Base.:^(z::Tv, n::AbstractFloat) where {Tv<:ADScaler}
    return Tv(z.v^n, z.i, z.∇*n*z.v^(n-1))
end




## SVar and SVar
function Base.:+(z::ADScaler{Nv, 1, Nv}, w::ADScaler{Nv, 1, Nv}) where {Nv}
    zi, wi = z.i[1], w.i[1]
    if zi == wi
        return ADScaler{Nv, 1, Nv}(z.v + w.v, z.i, z.∇ + w.∇)
    elseif wi < zi
        return ADScaler{Nv, 2, 2*Nv}(z.v + w.v, hcat(w.i, z.i), hcat(w.∇, z.∇))
    else
        return ADScaler{Nv, 2, 2*Nv}(z.v + w.v, hcat(z.i, w.i), hcat(z.∇, w.∇))
    end
end

function Base.:-(z::ADScaler{Nv, 1, Nv}, w::ADScaler{Nv, 1, Nv}) where {Nv}
    zi, wi = z.i[1], w.i[1]
    if zi == wi
        return ADScaler{Nv, 1, Nv}(z.v - w.v, z.i, z.∇ - w.∇)
    elseif wi < zi
        return ADScaler{Nv, 2, 2*Nv}(z.v - w.v, hcat(w.i, z.i), hcat(-w.∇, z.∇))
    else
        return ADScaler{Nv, 2, 2*Nv}(z.v - w.v, hcat(z.i, w.i), hcat(z.∇, -w.∇))
    end
end

function Base.:*(z::ADScaler{Nv, 1, Nv}, w::ADScaler{Nv, 1, Nv}) where {Nv}
    zi, wi = z.i[1], w.i[1]
    if zi == wi
        return ADScaler{Nv, 1, Nv}(z.v * w.v, z.i, w.v*z.∇ + z.v*w.∇)
    elseif wi < zi
        return ADScaler{Nv, 2, 2*Nv}(z.v * w.v, hcat(w.i, z.i), hcat(z.v*w.∇, w.v*z.∇))
    else
        return ADScaler{Nv, 2, 2*Nv}(z.v * w.v, hcat(z.i, w.i), hcat(w.v*z.∇, z.v*w.∇))
    end
end

function Base.:/(z::ADScaler{Nv, 1, Nv}, w::ADScaler{Nv, 1, Nv}) where {Nv}
    zv, wv = z.v, w.v
    zi, wi = z.i[1], w.i[1]
    if zi == wi
        return ADScaler{Nv, 1, Nv}(zv / wv, z.i,  (wv*z.∇ - zv*w.∇) / (wv*wv))
    elseif wi < zi
        return ADScaler{Nv, 2, 2*Nv}(zv / wv, hcat(w.i, z.i), hcat(-zv/(wv*wv)*w.∇, z.∇/wv))
    else
        return ADScaler{Nv, 2, 2*Nv}(zv / wv, hcat(z.i, w.i), hcat(z.∇/wv, -zv/(wv*wv)*w.∇))
    end
end

## SVar and DVar
function Base.:+(z::ADScaler{Nv, 1, Nv}, w::ADScaler{Nv, 2, L})::ADScaler{Nv, 2, L} where {Nv, L}
    g = zeros(SMat{Nv, 1, Nv})
    if z.i[1] == w.i[1]
        return ADScaler{Nv, 2, L}(z.v+w.v, w.i, hcat(z.∇, g) + w.∇)
    else
        return ADScaler{Nv, 2, L}(z.v+w.v, w.i, hcat(g, z.∇) + w.∇)
    end
end

function Base.:+(w::ADScaler{Nv, 2, L}, z::ADScaler{Nv, 1, Nv})::ADScaler{Nv, 2, L} where {Nv, L}
    g = zeros(SMat{Nv, 1, Nv})
    if z.i[1] == w.i[1]
        return ADScaler{Nv, 2, L}(z.v+w.v, w.i, hcat(z.∇, g) + w.∇)
    else
        return ADScaler{Nv, 2, L}(z.v+w.v, w.i, hcat(g, z.∇) + w.∇)
    end
end

function Base.:-(z::ADScaler{Nv, 1, Nv}, w::ADScaler{Nv, 2, L})::ADScaler{Nv, 2, L} where {Nv, L}
    g = zeros(SMat{Nv, 1, Nv})
    if z.i[1] == w.i[1]
        return ADScaler{Nv, 2, L}(z.v-w.v, w.i, hcat(z.∇, g) - w.∇)
    else
        return ADScaler{Nv, 2, L}(z.v-w.v, w.i, hcat(g, z.∇) - w.∇)
    end
end

function Base.:-(w::ADScaler{Nv, 2, L}, z::ADScaler{Nv, 1, Nv})::ADScaler{Nv, 2, L} where {Nv, L}
    g = zeros(SMat{Nv, 1, Nv})
    if z.i[1] == w.i[1]
        return ADScaler{Nv, 2, L}(w.v - z.v, w.i, w.∇ - hcat(z.∇, g))
    else
        return ADScaler{Nv, 2, L}(w.v - z.v, w.i, w.∇ - hcat(g, z.∇))
    end
end

function Base.:*(z::ADScaler{Nv, 1, Nv}, w::ADScaler{Nv, 2, L})::ADScaler{Nv, 2, L} where {Nv, L}
    g = zeros(SMat{Nv, 1, Nv})
    zv, wv = z.v, w.v
    if z.i[1] == w.i[1]
        return ADScaler{Nv, 2, L}(zv*wv, w.i, hcat(wv*z.∇, g) + zv*w.∇)
    else
        return ADScaler{Nv, 2, L}(zv*wv, w.i, hcat(g, wv*z.∇) + zv*w.∇)
    end
end

function Base.:*(w::ADScaler{Nv, 2, L}, z::ADScaler{Nv, 1, Nv})::ADScaler{Nv, 2, L} where {Nv, L}
    g = zeros(SMat{Nv, 1, Nv})
    zv, wv = z.v, w.v
    if z.i[1] == w.i[1]
        return ADScaler{Nv, 2, L}(zv*wv, w.i, hcat(wv*z.∇, g) + zv*w.∇)
    else
        return ADScaler{Nv, 2, L}(zv*wv, w.i, hcat(g, wv*z.∇) + zv*w.∇)
    end
end

function Base.:/(z::ADScaler{Nv, 1, Nv}, w::ADScaler{Nv, 2, L})::ADScaler{Nv, 2, L} where {Nv, L}
    g = zeros(SMat{Nv, 1, Nv})
    zv, wv = z.v, w.v
    if z.i[1] == w.i[1]
        return ADScaler{Nv, 2, L}(zv/wv, w.i, hcat(z.∇/wv, g) - zv/(wv*wv)*w.∇)
    else
        return ADScaler{Nv, 2, L}(zv/wv, w.i, hcat(g, z.∇/wv) - zv/(wv*wv)*w.∇)
    end
end

function Base.:/(w::ADScaler{Nv, 2, L}, z::ADScaler{Nv, 1, Nv})::ADScaler{Nv, 2, L} where {Nv, L}
    g = zeros(SMat{Nv, 1, Nv})
    zv, wv = z.v, w.v
    if z.i[1] == w.i[1]
        return ADScaler{Nv, 2, L}(wv/zv, w.i, hcat(-wv/(zv*zv)*z.∇, g) + w.∇/zv)
    else
        return ADScaler{Nv, 2, L}(wv/zv, w.i, hcat(g, -wv/(zv*zv)*z.∇) + w.∇/zv)
    end
end

## DVar and DVar
function Base.:+(z::ADScaler{Nv, 2, L}, w::ADScaler{Nv, 2, L}) where {Nv,L}
    # Limit to operations where z.i == w.i (z.i and w.i are always sorted)
    # That is, summation of two ADScalers belonging to the same connection
    return ADScaler{Nv, 2, L}(z.v+w.v, z.i, z.∇ + w.∇)
end

function Base.:-(z::ADScaler{Nv, 2, L}, w::ADScaler{Nv, 2, L}) where {Nv,L}
    # Limit to operations where z.i == w.i (z.i and w.i are always sorted)
    # That is, substraction of two ADScalers belonging to the same connection
    return ADScaler{Nv, 2, L}(z.v-w.v, z.i, z.∇ - w.∇)
end


## Function overload
Base.:(==)(z::Tv, w::Tv) where{Tv <: ADScaler}   = z.v == w.v
Base.:(==)(z::Tv, x::Real) where{Tv <: ADScaler}  = z.v == x
Base.:(==)(x::Real, z::Tv) where{Tv <: ADScaler}  = z.v == x

Base.isequal(z::Tv, w::Tv) where{Tv <: ADScaler} = isequal(z.v,w.v) && isequal(z.i, w.i) && isequal(z.∇, w.∇)
# Base.isequal(z::Tv, x::Real) where{Tv <: ADScaler}   = isequal(z.v, x) && isempty(z.∇)
# Base.isequal(x::Real, z::Tv) where{Tv <: ADScaler}   = isequal(z, x)

Base.isless(z::Tv,w::Tv) where{Tv <: ADScaler}  = z.v < w.v
Base.isless(z::Real,w::Tv) where{Tv <: ADScaler}  = z < w.v
Base.isless(z::Tv,w::Real) where{Tv <: ADScaler}  = z.v < w

Base.floor(z::Tv) where{Tv <: ADScaler}  = floor(z.v)
Base.ceil(z::Tv) where{Tv <: ADScaler}   = ceil(z.v)
Base.trunc(z::Tv) where{Tv <: ADScaler}  = trunc(z.v)
Base.round(z::Tv) where{Tv <: ADScaler}  = round(z.v)
Base.floor(::Type{T}, z::Tv) where {N, T<:Real, Tv <: ADScaler} = floor(T, z.v)
Base.trunc(::Type{T}, z::Tv) where {N, T<:Real, Tv <: ADScaler} = trunc(T, z.v)
Base.ceil( ::Type{T}, z::Tv) where {N, T<:Real, Tv <: ADScaler} = ceil(T, z.v)
Base.round(::Type{T}, z::Tv) where {N, T<:Real, Tv <: ADScaler} = round(T, z.v)

Base.abs(z::Tv) where{Tv <: ADScaler} = z ≥ 0 ? z : -z


end
