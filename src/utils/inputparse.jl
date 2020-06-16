module InputParse

function convert_float2vec(p::Union{Float64, Vector{Float64}}, nc::Int)
    return isa(p, Float64) ? p * ones(nc) : p
end

function parse_input_grid(input::Dict)::Dict
    grid_type = get(input, "grid", "Cartesian")
    if grid_type == "Cartesian"
        input["nc"] = input["nx"] * input["ny"] * input["nz"]
        v = ("dx", "dy", "dz")
        for p in v input[p] = convert_float2vec(input[p], input["nc"]) end
    end
    # Check required parameters
    v = ("d",)
    for p in v input[p] = convert_float2vec(input[p], input["nc"]) end
    return input
end

function parse_input_rock(input::Dict)::Dict
    # Check required parameters
    v = ("perm", "poro")
    for p in v input[p] = convert_float2vec(input[p], input["nc"]) end
    return input
end


function parse_input_fluid(input::Dict)::Dict
    state_type = input["fluid"]
    if state_type == "OW"
        v = ("po", "sw")
        for p in v input[p] = convert_float2vec(input[p], input["nc"]) end
    end
    return input
end

end
