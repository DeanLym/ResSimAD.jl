
const valid_well_modes = Dict([
    ("producer", ("bhp", "shut", "orat", "wrat", "lrat"))
    ("injector", ("bhp", "wrat", "shut"))
])

const valid_well_limits = keys(get_limit)

function parse_well_option(w, type::String, well_names::Set{String}, require_well_index::Bool)
    # Check required keywords
    wkeywords = keys(w)
    if !("name" ∈ wkeywords)
        error(LOGGER, "Missing well name")
    end
    vecs = ("perforation", "mode", "target")
    for pp in vecs
        if !(pp ∈ wkeywords)
            error(LOGGER, "Keyword \"$pp\" missing for well")
        end
    end

    # Parse well name
    check_keyword_type("name", w, (String,))
    name = w["name"]
    if name in well_names
        error(LOGGER, "Well name \"$name\" already exists")
    end
    push!(well_names, name)

    # Parse perforation
    check_keyword_type("perforation", w, (Vector{Tuple{Int, Int, Int}},))

    # Parse well control mode
    check_keyword_type("mode", w, (String,))
    mode = lowercase(w["mode"])
    valid_modes = valid_well_modes[type]
    if !(mode ∈ valid_modes)
        error(LOGGER, "Unsupported mode \"$mode\" for well \"$name\"")
    end
    w["mode"] = get_ctrl_mode[mode]

    # Parse well control target
    check_keyword_type("target", w, (Float64,))

    # Parse well radius
    if "wi" ∈ wkeywords
        check_keyword_type("wi", w, (Vector{Float64},))
        if length(w["wi"]) != length(w["perforation"])
            error(LOGGER, "Length of \"wi\" must be the same as \"perforation\" for well \"$name\"")
        end
        info(LOGGER, "Setting well index from input for well \"$name\"")
        if "radius" ∈ wkeywords
            warn(LOGGER, "Well radius will not be used for calculating well index for well \"$name\"")
        end
    elseif require_well_index
        error(LOGGER, "\"wi\" keyword must be set for well \"$name\"")
    else
        info(LOGGER, "Well index for well \"$name\" will be computed with radius")
        if !("radius" ∈ wkeywords)
            notice(LOGGER, "Keyword \"radius\" missing for well \"$name\", (default value 0.5 will be used)")
            w["radius"] = 0.5
        end
    end

    # Parse well limits
    if "limits" ∈ wkeywords
        check_keyword_type("limits", w, (Vector{Tuple{String, Float64}},))
        well_limits = copy(w["limits"])
        empty!(w["limits"])
        for (limit, v) in well_limits
            limit = lowercase(limit)
            if !(limit ∈ valid_well_limits)
                erorr(LOGGER, "Unsupported limit \"$limit\" for well \"$name\"")
            end
            if occursin(mode, limit)
                warn(LOGGER, "Limit \"$limit\" ignored under \"$mode\" control for well \"$name\"")
            else
                push!(w["limits"], (limit, v))
            end
        end
    else
        if type == "producer"
            w["limits"] = [("min_bhp", 14.7)]
        else
            w["limits"] = [("max_bhp", 1.0e5)]
        end
    end
end


function parse_facility(options, keywords, rock_opt)
    facility_opt = Dict()
    v = ("producers", "injectors")
    require_well_index = (rock_opt["trans_type"] == "trans")
    well_names = Set{String}()
    for p in v
        if p in keywords
            val = options[p]
            facility_opt["num_$p"] = length(val)
            facility_opt[p] = deepcopy(val)
            for w in facility_opt[p]
                parse_well_option(w, p[1:end-1], well_names, require_well_index)
            end
        else
            facility_opt["num_$p"] = 0
        end
    end
    return facility_opt
end
