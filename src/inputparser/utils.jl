function check_required_keyword(p, keywords)
    if !(p in keywords)
        error(LOGGER, "Keyword missing: $p")
    end
end


function check_keyword_type(p, options, types)
    type = typeof(options[p])
    if type <: Dict && Dict ∈ types
        return true
    elseif !(typeof(options[p]) ∈ types)
        error(LOGGER, "Type for keyword \"$p\" must be one of $types")
        return false
    end
    return true
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


function parse_value(x::String)::Tuple{Int, Float64}
    if '*' in x
        ind = findfirst('*', x)
        return parse(Int, x[1:ind-1]), parse(Float64, x[ind+1:end])
    else
        return 1, parse(Float64, x)
    end
end

function read_vector_keyword(fn::String, keyword::String, nc::Int)
    lines = filter!(x -> !startswith(x, "--"), readlines(fn))  # remove comments (line starting with -- )
    contents = strip(join(lines, " "))

    if !startswith(contents, keyword)
        error(LOGGER, "File $fn should start with $keyword ")
    end
    if !endswith(contents, "/")
        error(LOGGER, "File $fn should end with / ")
    end

    ind0 = findfirst(keyword, contents)[end] + 1 
    ind1 = findfirst("/", contents)[1] - 1

    values = split(contents[ind0+1:ind1])

    ret = zeros(Float64, nc)

    ct = 1
    for x in values
        n, v = parse_value(String(x))
        ret[ct:ct+n-1] .= v
        ct += n
    end

    return ret
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
        x = read_vector_keyword(options[p], uppercase(p), nc)
    else
        x = options[p]
    end
    x = check_dimension(x, nc, p)
    return x
end

function read_table_file(fn, keyword)
    contents = strip(open(f->read(f, String), fn))
    if !startswith(contents, keyword)
        error(LOGGER, "File $fn should start with $keyword ")
    end
    if !endswith(contents, "/")
        error(LOGGER, "File $fn should end with / ")
    end

    return DataFrame(CSV.File(IOBuffer(contents);
        delim=' ', comment="--", header=false, datarow=2, 
        footerskip=1, ignorerepeated=true, ignoreemptylines=true))
end