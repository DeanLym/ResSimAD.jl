abstract type AbstractPVT end


struct PVT <: AbstractPVT
    b::Interp
    μ::Interp
    table::DataFrame
end

function PVT(fn::String)
    # table = CSV.read(fn; delim=' ', comment="--", header=false, datarow=2, footerskip=1)
    # table = DataFrame(table)
    df = DataFrame!(CSV.File(fn; delim=' ', comment="--", header=false, datarow=2, footerskip=1))
    rename!(df, [:p, :b, :μ])
    # Add flat extrapolation
    insert!.(eachcol(df), 1, [-1.e30, df[1,:b], df[1,:μ])
    insert!.(eachcol(df), size(df)[1]+1, [1.e30, df[end,:b], df[end,:μ])
    #
    b = interpolate((df.p,), df.b, Gridded(Linear()))
    μ = interpolate((df.p,), df.μ, Gridded(Linear()))
    return PVT(b, μ, df)
end

function get_b(pvt::PVT, b::T, p::T)::T where {T <: AbstractVector}
    @. b = pvt.b(p)
end

function get_μ(pvt::PVT, μ::T, p::T)::T where {T <: AbstractVector}
    @. μ = pvt.μ(p)
end

struct PVTC <: AbstractPVT
    pref::Float64
    bref::Float64
    c::Float64
    μref::Float64
    cμ::Float64
end

function PVTC(param::Dict{String, Float64})
    vecs = ["pref", "bref",  "c", "μref", "cμ"]
    params = [param[v] for v in vecs]
    return PVTW(params...)
end

function PVTC(fn::String)
    # table = CSV.read(fn; delim=' ', comment="--", header=false, datarow=2, footerskip=1)
    # table = DataFrame(table)
    table = DataFrame!(CSV.File(fn; delim=' ', comment="--", header=false, datarow=2, footerskip=1))
    vecs = [:pref, :bref, :c, :μref, :cμ]
    rename!(table, vecs)
    params = [table[1, v] for v in vecs]
    return PVTC(params...)
end

function get_b(pvtc::PVTC, b::T, p::T)::T where {T <: AbstractVector}
    c, pref, bref = pvtc.c, pvtc.pref, pvtc.bref
    @. b = bref / (1 + c * (p - pref) + c*c*(p - pref)^2 / 2)
end

function get_μ(pvtc::PVTC, μ::T, p::T)::T where {T <: AbstractVector}
    cμ, pref, μref = pvtc.cμ, pvtc.pref, pvtc.μref
    @. μ = μref / (1 + cμ * (p - pref) + cμ * cμ * (p - pref)^2 / 2)
end
