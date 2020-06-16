abstract type AbstractPVT end


struct PVT <: AbstractPVT
    b::Interp
    μ::Interp
    table::DataFrame
end

function PVT(fn::String)
    table = CSV.read(fn; delim=' ', comment="--", header=false, datarow=2, footerskip=1)
    table = DataFrame(table)
    rename!(table, [:p, :b, :μ])
    b = interpolate((table.p,), table.b, Gridded(Linear()))
    μ = interpolate((table.p,), table.μ, Gridded(Linear()))
    return PVT(b, μ, table)
end

function get_b(pvt::PVT, b::Tensor, p::Tensor)::Tensor
    b .= pvt.b.(p)
end

function get_μ(pvt::PVT, μ::Tensor, p::Tensor)::Tensor
    μ .= pvt.μ.(p)
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
    table = CSV.read(fn; delim=' ', comment="--", header=false, datarow=2, footerskip=1)
    table = DataFrame(table)
    vecs = [:pref, :bref, :c, :μref, :cμ]
    rename!(table, vecs)
    params = [table[1, v] for v in vecs]
    return PVTC(params...)
end

function get_b(pvtc::PVTC, b::Tensor, p::Tensor)::Tensor
    c, pref, bref = pvtc.c, pvtc.pref, pvtc.bref
    x = c .* (p .- pref)
    b .= bref ./ (1 .+ x .+ x.^2 ./ 2)
end

function get_μ(pvtc::PVTC, μ::Tensor, p::Tensor)::Tensor
    cμ, pref, μref = pvtc.cμ, pvtc.pref, pvtc.μref
    x = cμ .* (p .- pref)
    μ .= μref ./ (1 .+ x .+ x.^2 ./ 2)
end
