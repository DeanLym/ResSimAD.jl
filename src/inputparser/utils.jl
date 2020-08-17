function check_required_keyword(p, keywords)
    if !(p in keywords)
        error(LOGGER, "Keyword missing: $p")
    end
end


function check_keyword_type(p, options, types)
    if !(typeof(options[p]) ∈ types)
        error(LOGGER, "Type for keyword \"$p\" must be one of $types")
    end
end


function check_dimension(
    p::Union{Float64,Vector{Float64}},
    n::Int,
    name::String;
    str = "all cells",
)
    if isa(p, Vector{Float64})
        if length(p) != n
            error(LOGGER, "Incorrect dimension for options $name")
        end
        return p
    else
        notice(LOGGER, "Setting constant $name $(round(p, digits=3)) for $str")
        return p * ones(n)
    end
end


function read_vector_keyword(fn::String, nc::Int)
    p = readdlm(fn; header=false, comments=true, comment_char='#')
    return reshape(p, nc)
end


function check_complementary_keywords(v, keywords)
    exist = map(x -> x in keywords, v)
    if !(any(exist))
        error(LOGGER, "Missing one of the keywords in $v")
    end
    v1 = v[findfirst(exist)]
    if sum(exist) > 1
        v0 = v[findall(exist)]
        v2 = v[2:end]
        warn(LOGGER, "Conflicting keywords: $v0")
        warn(LOGGER, "\"$v1\" will be used, ignoring $v2")
    end
    return v1
end


function parse_vector_keyword(p, options, nc)
    check_keyword_type(p, options, (Float64, Vector{Float64}, String))
    if isa(options[p], String)
        x = read_vector_keyword(options[p], nc)
    else
        x = options[p]
    end
    x = check_dimension(x, nc, p)
    return x
end
