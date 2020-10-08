
function parse_lsolver(options, keywords)
    lsolver_opt = Dict()
    p = "linear_solver_backend"
    if p ∈ keywords
        val = options[p]
        if val ∈ ("Julia", "DuneIstl")
            lsolver_opt["backend"] = options[p]
        else
            error(LOGGER, "Unknown linear solver backend $val")
        end
    else
        val = "DuneIstl"
        notice(LOGGER, "Keyword \"$p\" is missing (default value $val will be used)")
        lsolver_opt["backend"] = val
    end

    p = "linear_solver_type"
    if p ∈ keywords
        val = options[p]
        if lsolver_opt["backend"] == "Julia"
            supported = ("GMRES_ILU", "GMRES_CPR", "BICGSTAB_ILU", "BICGSTAB_CPR", "JULIA_BACKSLASH")
            if val ∈ supported
                lsolver_opt["type"] = options[p]
            else
                error(LOGGER, "Unsupported linear solver type $val by Julia backend. Supported linear solver types are $supported.")
            end
        else
            supported = ("GMRES_ILU", "BICGSTAB_ILU")
            if val ∈ supported
                lsolver_opt["type"] = options[p]
            else
                error(LOGGER, "Unsupported linear solver type $val by DuneIstl backend. Supported linear solver types are $supported.")
            end
        end
    else
        val = "BICGSTAB_ILU"
        notice(LOGGER, "Keyword \"$p\" is missing (default value $val will be used)")
        lsolver_opt["type"] = val
    end
    return lsolver_opt
end
