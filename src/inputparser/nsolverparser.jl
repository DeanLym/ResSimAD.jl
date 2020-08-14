
function parse_nsolver(options, keywords)
    nsolver_opt = Dict()
    v = ("max_newton_iter", "min_err")
    default_values = (10, 1.0e-6)
    for (p, val) in zip(v, default_values)
        if p ∈ keywords
            nsolver_opt[p] = options[p]
        else
            notice(LOGGER, "Keyword \"$p\" missing (default value $val will be used)")
            nsolver_opt[p] = val
        end
    end
    return nsolver_opt
end
