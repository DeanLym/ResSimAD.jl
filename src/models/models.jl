module Models

include("example1.jl")

function get_model(name::String)
    return @eval $(Symbol(name))()
end

end
