function parse_grid(options, keywords)
    # Detect grid type, only support Cartesian for now
    grid_type = "Cartesian"
    if grid_type == "Cartesian"
        info(LOGGER, "Grid type: Cartesian")
        return parse_cartesian_grid(options, keywords)
    end
end

function parse_cartesian_grid(options, keywords)
    grid_opt = Dict()
    grid_opt["type"] = "Cartesian"

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

    v = ("dx", "dy", "dz")
    for p in v check_required_keyword(p, keywords) end
    for p in v
        grid_opt[p] = parse_vector_keyword(p, options, nc)
    end
    v = ("tops", "d")
    p = check_complementary_keywords(v, keywords)
    if p == "tops"
        nx, ny, nz =  grid_opt["nx"], grid_opt["ny"], grid_opt["nz"]
        dz = grid_opt["dz"]
        tops = check_dimension(options[p], nx*ny, p; str="all top layer cells")
        grid_opt["d"] = compute_cell_depth(tops, dz, nx, ny, nz)
    else
        grid_opt["d"] = check_dimension(options[p], nc, p)
    end
    return grid_opt
end
