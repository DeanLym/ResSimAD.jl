module Models

using ..SimMaster: Sim

using DelimitedFiles

function get_example_data(name::String)
    return joinpath(@__DIR__, "data", name)
end

include("example1.jl")
include("example2.jl")

function get_model(name::String)
    return @eval $(Symbol(name))()
end


end
