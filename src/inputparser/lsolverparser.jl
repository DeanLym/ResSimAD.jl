
function parse_lsolver(options, keywords)
    lsolver_opt = Dict()
    p = "linear_solver_backend"
    if p ∈ keywords
        check_keyword_type(p, options, (String,))
        val = uppercase(options[p])
        if val ∈ ("JULIA", "DUNEISTL")
            lsolver_opt["backend"] = options[p]
        else
            error(LOGGER, "Unknown linear solver backend $val")
        end
    else
        val = "DUNEISTL"
        notice(LOGGER, "Keyword \"$p\" is missing (default value $val will be used)")
        lsolver_opt["backend"] = val
    end

    p = "linear_solver_type"
    if p ∈ keywords
        check_keyword_type(p, options, (String,))
        val = options[p]
        if lsolver_opt["backend"] == "JULIA"
            supported = ("GMRES_ILU", "GMRES_CPR", "BICGSTAB_ILU", "BICGSTAB_CPR", "JULIA_BACKSLASH")
            if val ∈ supported
                lsolver_opt["type"] = options[p]
            else
                error(LOGGER, "Unsupported linear solver type $val by JULIA backend. Supported linear solver types are $supported.")
            end
        else
            supported = ("GMRES_ILU", "BICGSTAB_ILU")
            if val ∈ supported
                lsolver_opt["type"] = options[p]
            else
                error(LOGGER, "Unsupported linear solver type $val by DUNEISTL backend. Supported linear solver types are $supported.")
            end
        end
    else
        val = "BICGSTAB_ILU"
        notice(LOGGER, "Keyword \"$p\" is missing (default value $val will be used)")
        lsolver_opt["type"] = val
    end
    return lsolver_opt
end
