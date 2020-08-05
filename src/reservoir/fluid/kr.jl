abstract type AbstractKROW end

struct SWOFTable <: AbstractKROW
    krw::Interp
    kro::Interp
    table::DataFrame
end

struct SWOFCorey <: AbstractKROW
    swi::Float64
    sor::Float64
    aw::Float64
    ao::Float64
    krw0::Float64
    kro0::Float64
    ds::Float64
end

function SWOFTable(fn::String)
    df = DataFrame(CSV.File(fn; delim=' ', comment="--", header=false, datarow=2, footerskip=1))
    rename!(df, [:sw, :krw, :kro, :pcw])
    # Add flat extrapolation
    insert!.(eachcol(df), 1, [-1.e30, df[1, :krw], df[1, :kro], df[1, :pcw]])
    insert!.(eachcol(df), size(df)[1]+1, [1.e30, df[end,:krw], df[end,:kro], df[end,:pcw]])
    #
    krw = interpolate((df.sw,), df.krw, Gridded(Linear()))
    kro = interpolate((reverse(1 .- df.sw),), reverse(df.kro), Gridded(Linear()))
    return SWOFTable(krw, kro, df)
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

function get_krw(krow::SWOFCorey, krw::ADVector, sw::ADVector)::ADVector
    krw0, ds, swl, swr, aw = krow.krw0, krow.ds, krow.swi, 1 - krow.sor, krow.aw
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

function get_kro(krow::SWOFCorey, kro::ADVector, so::ADVector)::ADVector
    kro0, ds, sor, sol, ao = krow.kro0, krow.ds, krow.sor, 1 - krow.swi, krow.ao
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

function get_krw(krow::SWOFTable, krw::ADVector, sw::ADVector)::ADVector
    krw .= krow.krw.(sw)
end

function get_kro(krow::SWOFTable, kro::ADVector, so::ADVector)::ADVector
    kro .= krow.kro.(so)
end
