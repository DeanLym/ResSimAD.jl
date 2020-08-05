module Models

include("example1.jl")

function get_model(name::String)
    return @eval $(Symbol(name))()
end

function get_table(name::String)
    return joinpath(@__DIR__, "tables", name)
end

end
