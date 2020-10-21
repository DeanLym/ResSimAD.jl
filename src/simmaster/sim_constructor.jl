const get_well_type =
    Dict{String,WellType}("producer" => PRODUCER, "injector" => INJECTOR)

function init_grid(grid_opt::Dict)
    info(LOGGER, "Initializing grid")
    grid_type = grid_opt["type"]
    if grid_type == "CARTESIAN"
        nx, ny, nz = grid_opt["nx"], grid_opt["ny"], grid_opt["nz"]
        grid = CartGrid(nx, ny, nz)
        if "dx" in keys(grid_opt)
            dx, dy, dz = grid_opt["dx"], grid_opt["dy"], grid_opt["dz"]
            set_cell_size(grid, dx, dy, dz)
        else
            set_cell_size(grid, grid_opt["v"])
        end
        set_cell_depth(grid, grid_opt["d"])
        info(LOGGER, "Grid dimension: $(grid.nx), $(grid.ny), $(grid.nz)")
        info(LOGGER, "Number of cells: $(grid.nc)")
        return grid
    else
        error(LOGGER, "Unsupported grid type $grid_type")
    end
end


function init_rock(rock_opt::Dict, nc::Int)
    info(LOGGER, "Initializing rock")
    if rock_opt["trans_type"] == "perm"
        rock = StandardRock(nc)
        set_perm(rock, rock_opt["permx"], rock_opt["permy"], rock_opt["permz"])
    else # rock_opt["trans_type"] == "trans"
        rock = TransRock(nc)
        set_trans(rock, rock_opt["tranx"], rock_opt["trany"], rock_opt["tranz"])
    end
    set_poro(rock, rock_opt["poro"])
    return rock
end

function init_fluid(fluid_opt::Dict, grid::AbstractGrid)
    info(LOGGER, "Initializing fluid")
    nc, nconn = grid.nc, grid.connlist.nconn

    fluid_type = fluid_opt["type"]
    if fluid_type == "OW"
        info(LOGGER, "Fluid system: $fluid_type")

        ρo, ρw = fluid_opt["ρo"], fluid_opt["ρw"]

        pvto_type, pvto_input = fluid_opt["PVTO_TYPE"], fluid_opt["PVTO"]
        if pvto_type == "PVDO"
            pvto = PVT(pvto_input)
        else
            pvto = PVTC(pvto_input)
        end

        pvtw = PVTC(fluid_opt["PVTW"])

        swof_type, swof_input = fluid_opt["SWOF_TYPE"], fluid_opt["SWOF"]

        if swof_type == "SWOF"
            swof = SWOFTable(swof_input)
        else
            swof = SWOFCorey(swof_input)
        end
        fluid = OWFluid(nc, nconn, ρo, ρw, pvto, pvtw, swof)
        # Set initial state
        sw = fluid_opt["sw"]
        po_type = fluid_opt["po_type"]
        if po_type == "po"
            po = fluid_opt["po"]
        else # po_type == "equil"
            # compute po from equil
            dref, pref = fluid_opt["equil"]
            Δd = mean(grid.dz) / 5.
            po = equil(fluid, dref, pref, grid.d, Δd)
        end
        set_fluid_tn(fluid, po, sw)
        update_primary_variable(fluid, po, sw)
        return fluid
    else
        error(LOGGER, "Unsupported fluid type $fluid_type")
    end
end


function init_well(T::WellType, w::Dict, nv::Int, grid::CartGrid, rock::AbstractRock)
    name = w["name"]
    perf = Vector{Int}()
    for indices in w["perforation"]
        push!(perf, get_grid_index(grid, indices...))
    end
    
    radius = w["radius"]
    well = StandardWell{T}(name, perf, radius, nv)
    # Set perforation depth
    well.d .= grid.d[perf]

    # Set control mode and target
    well.mode = w["mode"]
    well.target = w["target"]
    # Set limits
    if "limits" in keys(w)
        limits = w["limits"]
        for (k,v) in limits
            well.limits[get_limit[k]] = v
        end
    end
    compute_wi(well, grid, rock)
    return well
end

function init_well(T::WellType, w::Dict, nv::Int, grid::CartGrid )
    name = w["name"]
    perf = Vector{Int}()
    for indices in w["perforation"]
        push!(perf, get_grid_index(grid, indices...))
    end
    well = StandardWell{T}(name, perf, nv)
    # Set perforation depth
    well.d .= grid.d[perf]
    # Set control mode and target
    well.mode = w["mode"]
    well.target = w["target"]
    # Set limits
    if "limits" in keys(w)
        limits = w["limits"]
        for (k,v) in limits
            well.limits[get_limit[k]] = v
        end
    end
    # Set well index
    well.wi .= w["wi"]
    
    return well
end

function init_facility(facility_opt::Dict, nv::Int, grid::CartGrid, rock::StandardRock)
    info(LOGGER, "Initializing facility")
    v = ("producers", "injectors")
    facility = Dict{String, AbstractFacility}()
    for p in v
        T = get_well_type[p[1:end-1]]
        if facility_opt["num_$p"] > 0
            for w in facility_opt[p]
                facility[w["name"]] = init_well(T, w, nv, grid, rock)
            end
        end
    end

    return facility
end

function init_facility(facility_opt::Dict, nv::Int, grid::CartGrid, ::TransRock)
    info(LOGGER, "Initializing facility")
    v = ("producers", "injectors")
    facility = Dict{String, AbstractFacility}()
    for p in v
        T = get_well_type[p[1:end-1]]
        if facility_opt["num_$p"] > 0
            for w in facility_opt[p]
                facility[w["name"]] = init_well(T, w, nv, grid)
            end
        end
    end

    return facility
end


function init_scheduler(scheduler_opt::Dict)
    info(LOGGER, "Initializing scheduler")
    sch = Scheduler()
    sch.t_current = scheduler_opt["t_current"]
    sch.dt_max = scheduler_opt["dt_max"]
    sch.dt0 = scheduler_opt["dt0"]
    reset_dt(sch)
    set_time_step(sch, scheduler_opt["time_step"])
    return sch
end

function init_nsolver(nsolver_opt::Dict, lsolver_opt::Dict)
    info(LOGGER, "Initializing nonlinear solver")
    if lsolver_opt["backend"] == "Julia"
        nsolver = NRSolver()
    else
        nsolver = NRSolverDuneIstl()
    end
    nsolver.max_newton_iter = nsolver_opt["max_newton_iter"]
    nsolver.min_err =  nsolver_opt["min_err"]
    return nsolver
end

function init_lsolver(lsolver_opt::Dict, grid::AbstractGrid, fluid::AbstractFluid)
    info(LOGGER, "Initializing linear solver")
    solver_type = lsolver_opt["type"]
    solver_backend = lsolver_opt["backend"]
    if solver_backend == "Julia"
        if solver_type == "GMRES_ILU"
            lsolver = GMRES_ILU_Solver()
        elseif solver_type == "BICGSTAB_ILU"
            lsolver = BICGSTAB_ILU_Solver()
        elseif solver_type == "GMRES_CPR"
            lsolver = GMRES_CPR_Solver(grid.nc, grid.neighbors)
        elseif solver_type == "BICGSTAB_CPR"
            lsolver = BICGSTAB_CPR_Solver(grid.nc, grid.neighbors)
        else
            lsolver = Julia_BackSlash_Solver()
        end
    else
        lsolver = init_duneistl_solver(solver_type, grid.nc, fluid)
    end
end

function init_duneistl_solver(solver_type::String, nc::Int, ::OWFluid)
    solver = DuneIstlSolver{Float64, Int32(2)}(nc)
    if solver_type == "GMRES_ILU"
        set_solver_type(solver, "RestartedGMRes")
        set_preconditioner_type(solver, "ILU")
    else # "BICGSTAB_ILU"
        set_solver_type(solver, "BiCGSTAB")
        set_preconditioner_type(solver, "ILU")
    end
    return DuneIstlSolverWrapper(solver, Int[])
end

function init_duneistl_solver(solver_type::String, nc::Int, ::SPFluid)
    lsolver = DuneIstlSolver{Float64, Int32(1)}(nc)
    if solver_type == "GMRES_ILU"
        set_solver_type(lsolver, "RestartedGMRes")
        set_preconditioner_type(lsolver, "ILU")
    else # "BICGSTAB_ILU"
        set_solver_type(lsolver, "BiCGSTAB")
        set_preconditioner_type(lsolver, "ILU")
    end
    return lsolver
end


function Sim(options::Dict)
    info(LOGGER, "Creating simulation model")
    parsed_options = parse_input(options)

    grid = init_grid(parsed_options["grid_opt"])

    nc = grid.nc
    rock = init_rock(parsed_options["rock_opt"], nc)
    
    construct_connlist(grid, rock)
    sort_conn(grid.connlist)
    info(LOGGER, "Number of connections: $(grid.connlist.nconn)")

    nconn = grid.connlist.nconn
    fluid = init_fluid(parsed_options["fluid_opt"], grid)

    reservoir = StandardReservoir(grid, rock, fluid)

    nv = fluid.nv
    facility = init_facility(parsed_options["facility_opt"], nv, grid, rock)

    scheduler = init_scheduler(parsed_options["scheduler_opt"])

    nsolver = init_nsolver(parsed_options["nsolver_opt"], parsed_options["lsolver_opt"])

    lsolver = init_lsolver(parsed_options["lsolver_opt"], grid, fluid)

    nsolver.lsolver = lsolver

    update_phases(fluid, grid.connlist)

    compute_residual(fluid, grid, rock, facility, scheduler.dt)

    initialize_nsolver(nsolver, grid, fluid)

    return Sim(reservoir, facility, scheduler, nsolver)
end
