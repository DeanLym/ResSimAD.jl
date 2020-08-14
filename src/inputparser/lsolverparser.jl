
function parse_lsolver(options, keywords)
    lsolver_opt = Dict()
    p = "linear_solver"
    if p ∈ keywords
        val = options[p]
        if val ∈ ("GMRES_ILU0", "Julia_BackSlash")
            lsolver_opt["type"] = options[p]
        else
            error(LOGGER, "Unknown linear solver type $val")
        end
    else
        val = "GMRES_ILU0"
        notice(LOGGER, "Keyword \"$p\" is missing (default value $val will be used)")
        lsolver_opt["type"] = val
    end
    if lsolver_opt["type"] == "GMRES_ILU0"
        p = "τ"
        if p ∈ keywords
            val = options[p]
            if val < 0.0 || val > 1.0
                error(LOGGER, "Value for \"$p\" must be within (0.0, 1.0)")
            end
            lsolver_opt[p] = val
        else
            val = 0.1
            notice(LOGGER, "Keyword \"$p\" for GMRES_ILU0 is missing (default value $val will be used)")
            lsolver_opt[p] = val
        end
    end
    return lsolver_opt
end
