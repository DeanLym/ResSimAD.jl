module Fluid

using ..AutoDiff: Tensor


# abstract type AbstractSWOF end
# abstract type AbstractPVTO end
# abstract type AbstractPVTW end

# function compute_bo(pvdo::PVDO, bo::Tensor, po::Tensor)::Tensor
#     bo .= pvdo.bo.(po)
# end
#
# function compute_μo(pvdo::PVDO, μo::Tensor, po::Tensor)::Tensor
#     μo .= pvdo.μo.(po)
# end
#
# function compute_bw(pvtw::PVTW, bw::Tensor, pw::Tensor)::Tensor
#     cw, pw_ref, bw_ref = pvtw.cw, pvtw.pw_ref, pvtw.bw_ref
#     x = cw .* (pw .- pw_ref)
#     bw .= bw_ref ./ (1 .+ x .+ x.^2 ./ 2)
# end
#
# function compute_μw(pvtw::PVTW, μw::Tensor, pw::Tensor)::Tensor
#     cμ, pw_ref, μw_ref = pvtw.cμ, pvtw.pw_ref, pvtw.μw_ref
#     x = -cμ .* (pw .- pw_ref)
#     μw .= μw_ref ./ (1 .+ x .+ x.^2 ./ 2)
# end



#
# function PVDO(fn::String)
#     table = CSV.read(fn; delim=' ', comment="--", header=false, datarow=2, footerskip=1)
#     table = DataFrame(table)
#     rename!(table, [:po, :bo, :μo])
#     bo = interpolate((table.po,), table.bo, Gridded(Linear()))
#     μo = interpolate((table.po,), table.μo, Gridded(Linear()))
#     return PVDO(bo, μo, table)
# end
#
# function PVTW(param::Dict{String, Float64})
#     vecs = ["pw_ref", "bw_ref",  "cw", "μw_ref", "cμ"]
#     params = [param[v] for v in vecs]
#     return PVTW(params...)
# end
#
# function PVTW(fn::String)
#     table = CSV.read(fn; delim=' ', comment="--", header=false, datarow=2, footerskip=1)
#     table = DataFrame(table)
#     vecs = [:pw_ref, :bw_ref, :cw, :μw_ref, :cμ]
#     rename!(table, vecs)
#     params = [table[1, v] for v in vecs]
#     return PVTW(params...)
# end
#
# function PVTW(table::DataFrame)
#     vecs = [:pw_ref, :bw_ref, :cw, :μw_ref, :cμ]
#     params = [table[1, v] for v in vecs]
#     return PVTW(params...)
# end
#
# function PVDO(table::DataFrame)
#     bo = interpolate((table.po,), table.bo, Gridded(Linear()))
#     μo = interpolate((table.po,), table.μo, Gridded(Linear()))
#     return PVDO(bo, μo, table)
# end

#
# struct PVDO <: AbstractPVT
#     bo::Interp
#     μo::Interp
#     table::DataFrame
# end
#
# struct PVTW <: AbstractPVT
#     pw_ref::Float64
#     bw_ref::Float64
#     cw::Float64
#     μw_ref::Float64
#     cμ::Float64
# end


end
