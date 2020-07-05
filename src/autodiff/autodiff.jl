module AutoDiff

using StaticArrays
import Base
const SVec = SVector{Nv, Float64} where {Nv}
const SMat = SMatrix{Nv, Nc, Float64, L} where {Nv, Nc, L}

struct ADScaler{Nv, Nc, L}
    v::Float64
    i::SVector{Nc, Int}
    ∇::SMat{Nv, Nc, L}
end

const ADVector = Vector{T} where T <: ADScaler

## Constructors
function advector(v::Vector{Float64}, g::SMat{Nv, 1, L}) where{Nv, L}
    return [ADScaler{Nv, 1, L}(v[i], SA[i], g) for i in eachindex(v)]
end

function advector(v::Vector{Float64}, iv::Int, nv::Int)
    g = SMat{nv, 1, nv}([if i==iv 1.0 else 0.0 end for i=1:nv])
    return new_param(v, g)
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

function Base.:^(z::Tv, n::T) where {Tv<:ADScaler, T<:Real}
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


end
