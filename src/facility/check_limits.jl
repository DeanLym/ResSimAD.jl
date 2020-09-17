## check limits for producer

function check_limit(well::StandardWell{T}, limit::Limit, x::Float64, y::Float64, new_limit::Limit, new_mode::CtrlMode, unit::String) where {T}
    if x < y
        well.limits[new_limit] = well.target
        well.mode = new_mode
        well.target = well.limits[limit]
        warn(LOGGER, "Well \"$(well.name)\" reached $limit limit and switched to $new_mode at $(well.target) $unit")
        return false
    end
    return true
end

function check_limits(well::StandardWell{PRODUCER})::Bool
    mode = well.mode
    if mode == CBHP
        if MAX_ORAT in keys(well.limits)
            limit_honored = check_limit(well, MAX_ORAT, well.limits[MAX_ORAT], sum(value(well.qo)), MIN_BHP, CORAT, "STB/Day") 
            if !limit_honored
                return false
            end
        end
        if MAX_LRAT in keys(well.limits)
            limit_honored = check_limit(well, MAX_LRAT, well.limits[MAX_LRAT], sum(value(well.qo + well.qw)), MIN_BHP, CLRAT, "STB/Day") 
            if !limit_honored
                return false
            end
        end
    elseif mode == CORAT
        if MIN_BHP in keys(well.limits)
            limit_honored = check_limit(well, MIN_BHP, mean(value(well.bhp)), well.limits[MIN_BHP], MAX_ORAT, CBHP, "psi") 
            if !limit_honored
                return false
            end
        end
        if MAX_LRAT in keys(well.limits)
            imit_honored = check_limit(well, MAX_LRAT, well.limits[MAX_LRAT], sum(value(well.qo + well.qw)), MAX_ORAT, CLRAT, "STB/Day") 
            if !limit_honored
                return false
            end
        end
    elseif mode == CLRAT
        if MIN_BHP in keys(well.limits)
            limit_honored = check_limit(well, MIN_BHP, mean(value(well.bhp)), well.limits[MIN_BHP], MAX_LRAT, CBHP, "psi") 
            if !limit_honored
                return false
            end
        end
        if MAX_ORAT in keys(well.limits)
            limit_honored = check_limit(well, MAX_ORAT, well.limits[MAX_ORAT], sum(value(well.qo)), MAX_LRAT, CORAT, "STB/Day") 
            if !limit_honored
                return false
            end
        end
    end
    return true
end

## check limits for injector
function check_limits(well::StandardWell{INJECTOR})::Bool
    mode = well.mode
    if mode == CBHP
        if MAX_WRAT in keys(well.limits)
            limit_honored = check_limit(well, MAX_WRAT, sum(value(well.qw)), well.limits[MAX_WRAT], MAX_BHP, CWRAT, "STB/Day") 
            if !limit_honored
                return false
            end
        end
    elseif mode == CWRAT
        if MAX_BHP in keys(well.limits)
            limit_honored = check_limit(well, MAX_BHP, well.limits[MAX_BHP], mean(value(well.bhp)), MAX_WRAT, CBHP, "psi") 
            if !limit_honored
                return false
            end
        end
    end
    return true
end

