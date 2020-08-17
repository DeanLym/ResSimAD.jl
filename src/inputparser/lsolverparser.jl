
function parse_lsolver(options, keywords)
    lsolver_opt = Dict()
    p = "linear_solver"
    if p ∈ keywords
        val = options[p]
        if val ∈ ("GMRES_ILU0", "GMRES_CPR", "Julia_BackSlash")
            lsolver_opt["type"] = options[p]
        else
            error(LOGGER, "Unknown linear solver type $val")
        end
    else
        val = "GMRES_ILU0"
        notice(LOGGER, "Keyword \"$p\" is missing (default value $val will be used)")
        lsolver_opt["type"] = val
    end
    return lsolver_opt
end
