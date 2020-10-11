
function parse_grid(options, keywords)
    # Detect grid type, only support Cartesian for now
    grid_type = uppercase(get(options, "grid_type", "CARTESIAN"))
    if grid_type == "CARTESIAN"
        info(LOGGER, "Grid type: Cartesian")
        return parse_cartesian_grid(options, keywords)
    end
end

function parse_cartesian_grid(options, keywords)
    grid_opt = Dict()
    grid_opt["type"] = "CARTESIAN"
    # Check required keywords for setting grid dimension
    v = ("nx", "ny", "nz")
    for p in v check_required_keyword(p, keywords) end
    for p in v check_keyword_type(p, options, (Int64,)) end
    # Set grid dimension
    grid_opt["nc"] = 1
    for p in v
        grid_opt[p] = options[p]
        grid_opt["nc"] *= grid_opt[p]
    end
    nc = grid_opt["nc"]

    if "nc" ∈ keywords
        warn(LOGGER, "Ignoring keyword \"nc\" for Cartesian grid.")
    end

    # Parse dx, dy, dz
    v = ("dx", "dy", "dz")
    for p in v
        if p ∈ keywords
            grid_opt[p] = parse_vector_keyword(p, options, nc)
        end
    end
    
    # Check required keywords for calculating volume
    p1 = "v"
    if p1 ∈ keywords
        grid_opt[p1] = parse_vector_keyword(p1, options, nc)
    else
        info(LOGGER, "Keyword \"v\" not found, calculating volume from dx, dy, dz")
        for p in v # check dx, dy, dz are all specified
            check_required_keyword(p, keywords)
        end
        grid_opt["v"] = grid_opt["dx"] .* grid_opt["dy"] .* grid_opt["dz"]
    end

    # Check required keywords for calculating depth
    v = ("d", "tops")
    p = check_complementary_keywords(v, keywords)
    if p == "d"
        grid_opt["d"] = parse_vector_keyword(p, options, nc)
    else
        nx, ny, nz, dz = grid_opt["nx"], grid_opt["ny"], grid_opt["nz"], grid_opt["dz"]
        tops = check_dimension(options[p], nx*ny, p; str="all top layer cells")
        grid_opt["d"] = compute_cell_depth(tops, dz, nx, ny, nz)
    end

    return grid_opt
end
