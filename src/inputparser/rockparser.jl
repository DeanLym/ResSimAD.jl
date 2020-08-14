
function parse_rock(options, keywords, nc)
    rock_opt = Dict()
    p = "poro"
    check_required_keyword(p, keywords)
    rock_opt[p] = parse_vector_keyword(p, options, nc)
    # Parse permeability inputs
    p1 = "perm"
    v2 = ("permx", "permy", "permz")
    if p1 ∈ keywords
        exists = map(x-> x ∈ keywords, v2)
        if any(exists)
            ind = findall(exists)
            warn(LOGGER, "Conflicting keywords: $p1 and $(v2[ind])")
            warn(LOGGER, "\"$p1\" will be used, ignoring $(v2[ind])")
        end
        perm = parse_vector_keyword(p1, options, nc)
        for p in v2
            rock_opt[p] = copy(perm)
        end
    else
        for p in v2
            check_required_keyword(p, keywords)
            rock_opt[p] = parse_vector_keyword(p, options, nc)
        end
    end
    p = "multpermz"
    if p in keywords
        check_keyword_type(p, options, (Float64,))
        multz = options[p]
        if multz ≤ 0
            error(LOGGER, "Negative value for $p")
        end
        rock_opt["permz"] .*= multz
    end
    return rock_opt
end
