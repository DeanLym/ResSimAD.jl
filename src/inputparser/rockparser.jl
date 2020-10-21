
function parse_rock(options, keywords, nc, grid_opt)
    rock_opt = Dict()
    # Parse poro
    p = "poro"
    check_required_keyword(p, keywords)
    rock_opt[p] = parse_vector_keyword(p, options, nc)
    # Check required keywords for calculating transmissibility
    v = ("tranx", "trany", "tranz")
    exists = map(x-> x ∈ keywords, v)
    if any(exists)
        info(LOGGER, "Setting transmissibility from tranx, trany, tranz")
        rock_opt["trans_type"] = "trans" # set trans from input
        for p in v
            check_required_keyword(p, keywords)
            rock_opt[p] = parse_vector_keyword(p, options, nc)
        end
        for p in ("perm", "permx", "permy", "permz", "multpermz")
            if p ∈ keywords
                warn(LOGGER, "Conflicting keywords: $p and $v")
                warn(LOGGER, "Ignoring keyword $p")
            end
        end
    else
        rock_opt["trans_type"] = "perm" # compute trans from perm
        v2 = ("permx", "permy", "permz")
        exists = map(x-> x ∈ keywords, v2)
        if any(exists)
            info(LOGGER, "Setting permeability from permx, permy, permz")
            for p in v2
                check_required_keyword(p, keywords)
                rock_opt[p] = parse_vector_keyword(p, options, nc)
            end
            p = "perm"
            if p ∈ keywords
                warn(LOGGER, "Conflicting keywords: $p and $v2")
                warn(LOGGER, "Ignoring keyword $p")
            end
        elseif "perm"  ∈ keywords
            info(LOGGER, "Setting permeability from perm")
            perm = parse_vector_keyword("perm", options, nc)
            for p in v2
                rock_opt[p] = copy(perm)
            end
        else
            error(LOGGER, "Missing keywords for calculating transmissibility, need one of (perm, (permx, permy, permz), (tranx, trany, tranz))")
        end
        p = "multpermz"
        if p in keywords
            check_keyword_type(p, options, (Float64,))
            info(LOGGER, "Applying permz multiplier")
            multz = options[p]
            if multz ≤ 0
                error(LOGGER, "Negative value for $p")
            end
            rock_opt["permz"] .*= multz
        end
        # Check if dx, dy, dz are available in grid_opt
        for p in ("dx", "dy", "dz")
            if !(p∈keys(grid_opt))
                error(LOGGER, "When calculating transmissibility from permx, permy, permz, keyword \"$p\" is required")
            end
        end
    end
    return rock_opt
end
