
function parse_scheduler(options, keywords)
    scheduler_opt = Dict()
    scheduler_opt["t_current"] = 0.0
    v = ("dt0", "dt_max", "t_end")
    default_values = (0.1, 50., 100.)
    for (p, val) in zip(v, default_values)
        if p ∈ keywords
            scheduler_opt[p] = options[p]
        else
            notice(LOGGER, "Keyword \"$p\" missing (default value $val will be used)")
            scheduler_opt[p] = val
        end
    end
    p = "time_step"
    if p ∈ keywords
        check_keyword_type(p, options, (Vector{Float64},))
        scheduler_opt[p] = options[p]
    else
        val = [scheduler_opt["t_end"]]
        notice(LOGGER, "Keyword \"$p\" missing (default value $val will be used)")
        scheduler_opt[p] = val
    end
    return scheduler_opt
end
