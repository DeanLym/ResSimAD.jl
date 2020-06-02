module Fluid

using ..AutoDiff: Tensor

using DataFrames
using CSV
using Interpolations
using Interpolations:GriddedInterpolation

const Interp = GriddedInterpolation{Float64,1,Float64,Gridded{Linear},Tuple{Array{Float64,1}}}

abstract type AbstractSWOF end
abstract type AbstractPVTO end
abstract type AbstractPVTW end

struct PVDO <: AbstractPVTO
    bo::Interp
    μo::Interp
    table::DataFrame
end

struct PVTW <: AbstractPVTW
    pw_ref::Float64
    bw_ref::Float64
    cw::Float64
    μw_ref::Float64
    cμ::Float64
end

struct SWOFTable <: AbstractSWOF
    krw::Interp
    kro::Interp
    table::DataFrame
end

struct SWOFCorey <: AbstractSWOF
    swi::Float64
    sor::Float64
    aw::Float64
    ao::Float64
    krw0::Float64
    kro0::Float64
    ds::Float64
end

function PVDO(fn::String)
    table = CSV.read(fn; delim=' ', comment="--", header=false, datarow=2, footerskip=1)
    table = DataFrame(table)
    rename!(table, [:po, :bo, :μo])
    bo = interpolate((table.po,), table.bo, Gridded(Linear()))
    μo = interpolate((table.po,), table.μo, Gridded(Linear()))
    return PVDO(bo, μo, table)
end

function PVTW(param::Dict{String, Float64})
    vecs = ["pw_ref", "bw_ref",  "cw", "μw_ref", "cμ"]
    params = [param[v] for v in vecs]
    return PVTW(params...)
end

function PVTW(fn::String)
    table = CSV.read(fn; delim=' ', comment="--", header=false, datarow=2, footerskip=1)
    table = DataFrame(table)
    vecs = [:pw_ref, :bw_ref, :cw, :μw_ref, :cμ]
    rename!(table, vecs)
    params = [table[1, v] for v in vecs]
    return PVTW(params...)
end

function PVTW(table::DataFrame)
    vecs = [:pw_ref, :bw_ref, :cw, :μw_ref, :cμ]
    params = [table[1, v] for v in vecs]
    return PVTW(params...)
end

function PVDO(table::DataFrame)
    bo = interpolate((table.po,), table.bo, Gridded(Linear()))
    μo = interpolate((table.po,), table.μo, Gridded(Linear()))
    return PVDO(bo, μo, table)
end

function SWOFTable(fn::String)
    table = CSV.read(fn; delim=' ', comment="--", header=false, datarow=2, footerskip=1)
    table = DataFrame(table)
    rename!(table, [:sw, :krw, :kro, :pcw])
    krw = interpolate((table.sw,), table.krw, Gridded(Linear()))
    kro = interpolate((reverse(1 .- table.sw),), reverse(table.kro), Gridded(Linear()))
    return SWOFTable(krw, kro, table)
end


function SWOFTable(table::DataFrame)
    krw = interpolate((table.sw[:],), table.krw[:], Gridded(Linear()))
    kro = interpolate((reverse(1 .- table.sw[:]),), reverse(table.kro[:]), Gridded(Linear()))
    return SWOFTable(krw, kro, table)
end


function SWOFCorey(param::Dict{String, Float64})
    vecs = ["swi", "sor", "aw", "ao", "krw0", "kro0"]
    params = [param[v] for v in vecs]
    ds = 1 - param["swi"] - param["sor"]
    return SWOFCorey(params..., ds)
end


function compute_krw(swof::SWOFCorey, krw::Tensor, sw::Tensor)::Tensor
    krw0, ds, swl, swr, aw = swof.krw0, swof.ds, swof.swi, 1 - swof.sor, swof.aw
    for i = 1:length(sw)
        if sw[i] <= swl
            krw[i] = 0.0 * sw[i]  #To keep track of gradients
        elseif sw[i] >= swr
            krw[i] = 0.0 * sw[i] + krw0
        else
            krw[i] = krw0 * ((sw[i] - swi) / ds)^aw
        end
    end
    return krw
end

function compute_kro(swof::SWOFCorey, kro::Tensor, so::Tensor)::Tensor
    kro0, ds, sor, sol, ao = swof.kro0, swof.ds, swof.sor, 1 - swof.swi, swof.ao
    for i = 1:length(sw)
        if so[i] <= sor
            kro[i] = 0.0 * so[i]  #To keep track of gradients
        elseif so[i] >= sol
            kro[i] = 0.0 * so[i] + kro0
        else
            kro[i] = kro0 * ((so[i] - sor) / ds)^ao
        end
    end
    return kro
end

function compute_krw(swof::SWOFTable, krw::Tensor, sw::Tensor)::Tensor
    krw .= swof.krw.(sw)
end

function compute_kro(swof::SWOFTable, kro::Tensor, so::Tensor)::Tensor
    kro .= swof.kro.(so)
end

function compute_bo(pvdo::PVDO, bo::Tensor, po::Tensor)::Tensor
    bo .= pvdo.bo.(po)
end

function compute_μo(pvdo::PVDO, μo::Tensor, po::Tensor)::Tensor
    μo .= pvdo.μo.(po)
end

function compute_bw(pvtw::PVTW, bw::Tensor, pw::Tensor)::Tensor
    cw, pw_ref, bw_ref = pvtw.cw, pvtw.pw_ref, pvtw.bw_ref
    x = cw .* (pw .- pw_ref)
    bw .= bw_ref ./ (1 .+ x .+ x.^2 ./ 2)
end

function compute_μw(pvtw::PVTW, μw::Tensor, pw::Tensor)::Tensor
    cμ, pw_ref, μw_ref = pvtw.cμ, pvtw.pw_ref, pvtw.μw_ref
    x = -cμ .* (pw .- pw_ref)
    μw .= μw_ref ./ (1 .+ x .+ x.^2 ./ 2)
end

end
