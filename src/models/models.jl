module Models

using ..SimMaster: Sim

using DelimitedFiles

function get_example_data(name::String)
    return joinpath(@__DIR__, "data", name)
end

include("example1.jl")
include("example2.jl")
include("example3.jl")
include("example4.jl")


"""
    get_model(example_name::String)

Get example model `example_name`.

# Examples
```jldoctest
julia> using ResSimAD: get_model

julia> sim, options = get_model("example1");

```

"""
function get_model(example_name::String)
    return @eval $(Symbol(example_name))()
end


end
