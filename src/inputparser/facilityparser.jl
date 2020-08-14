
const well_modes_dict = Dict([
    ("producer", ("bhp", "shut", "orat", "wrat", "lrat"))
    ("injector", ("bhp", "wrat"))
])

function parse_well_option(w, type::String, well_names::Set{String})
    modes = well_modes_dict[type]
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
    check_keyword_type("name", w, (String,))
    name = w["name"]
    if name in well_names
        error(LOGGER, "Well name \"$name\" already exists")
    end
    push!(well_names, name)

    check_keyword_type("perforation", w, (Vector{Tuple{Int, Int, Int}},))
    check_keyword_type("mode", w, (String,))
    mode = lowercase(w["mode"])
    if !(mode ∈ modes)
        erorr(LOGGER, "Unsupported mode \"$mode\" for well \"$name\"")
    end
    w["mode"] = get_ctrl_mode[mode]
    check_keyword_type("target", w, (Float64,))

    if !("radius" ∈ wkeywords)
        notice(LOGGER, "Keyword \"radius\" missing for well \"$name\", (default value 0.5 will be used)")
        w["radius"] = 0.5
    end
end


function parse_facility(options, keywords)
    facility_opt = Dict()
    v = ("producers", "injectors")

    well_names = Set{String}()
    for p in v
        if p in keywords
            val = options[p]
            facility_opt["num_$p"] = length(val)
            facility_opt[p] = deepcopy(val)
            for w in facility_opt[p]
                parse_well_option(w, p[1:end-1], well_names)
            end
        else
            facility_opt["num_$p"] = 0
        end
    end
    return facility_opt
end
