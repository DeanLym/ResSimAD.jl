abstract type AbstractPVT end


struct PVT <: AbstractPVT
    b
    bμ
    table::DataFrame
end

function PVT(df::DataFrame)
    b = interpolate((df.p,), df.b, Gridded(Linear()))
    bμ = interpolate((df.p,), df.b .* df.μ , Gridded(Linear()))
    return PVT(b, bμ, df)
end

function get_b(pvt::PVT, b::T, p::T)::T where {T <: AbstractVector}
    @. b = pvt.b(p)
end

function get_μ(pvt::PVT, μ::T, p::T, b::T)::T where {T <: AbstractVector}
    @. μ = pvt.bμ(p) / b
end

struct PVTC <: AbstractPVT
    pref::Float64
    bref::Float64
    c::Float64
    μref::Float64
    cμ::Float64
end

function PVTC(param::Dict{String, Float64})
    vecs = ("pref", "bref",  "c", "μref", "cμ")
    params = [param[v] for v in vecs]
    return PVTC(params...)
end

function PVTC(df::DataFrame)
    vecs = (:pref, :bref, :c, :μref, :cμ)
    params = [df[1, v] for v in vecs]
    return PVTC(params...)
end

function get_b(pvtc::PVTC, b::T, p::T)::T where {T <: AbstractVector}
    c, pref, bref = pvtc.c, pvtc.pref, pvtc.bref
    @. b = bref / (1 + c * (p - pref) + c^2*(p - pref)^2 / 2)
end

function get_μ(pvtc::PVTC, μ::T, p::T, b::T)::T where {T <: AbstractVector}
    c, cμ, pref, μref, bref = pvtc.c, pvtc.cμ, pvtc.pref, pvtc.μref, pvtc.bref
    @. μ = bref*μref / (1 + (c-cμ) * (p - pref) + (c-cμ)^2 * (p - pref)^2 / 2) / b
end

