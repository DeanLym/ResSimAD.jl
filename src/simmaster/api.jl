
get_well_type =
    Dict{String,WellType}("producer" => PRODUCER, "injector" => INJECTOR)


function init_grid(input::Dict)
    input = parse_input_grid(input)
    nx, ny, nz = input["nx"], input["ny"], input["nz"]
    grid = CartGrid(nx, ny, nz)
    dx, dy, dz = input["dx"], input["dy"], input["dz"]
    set_cell_size(grid, dx, dy, dz)
    set_cell_depth(grid, input["d"])
    return grid, input
end

function init_rock(input::Dict)
    input = parse_input_rock(input)
    rock = StandardRock(input["nc"])
    set_perm(rock, input["perm"])
    set_poro(rock, input["poro"])
    return rock, input
end

function init_fluid(input::Dict, grid::AbstractGrid)
    input = parse_input_fluid(input)
    fluid = input["fluid"]
    if fluid == "OW"
        return init_owfluid(input, grid)
    elseif fluid ∈ ["OIL", "WATER", "GAS"]
        return init_spfluid(input, grid)
    end
end

function init_spfluid(input::Dict, grid::AbstractGrid)
    if input["PVT"] == "PVDO.DAT"
        pvt = PVT(input["PVT"])
    elseif input["PVT"] == "PVCDO.DAT"
        pvt = PVTC(input["PVT"])
    end
    ρ = get(input, "ρ", 49.1)
    nc, nconn = grid.nc, grid.connlist.nconn
    fluid = SPFluid(nc, nconn, ρ, pvt; name=input["fluid"])
    return fluid, input
end

function init_owfluid(input::Dict, grid::AbstractGrid)
    # Construct Property Tables
    if "SWOF" in keys(input)
        krow = SWOFTable(input["SWOF"])
    elseif "SWOFCorey" in keys(input)
        krow = SWOFCorey(input["SWOFCorey"])
    end
    pvto = PVT(input["PVDO"])
    pvtw = PVTC(input["PVTW"])
    # Construct State
    ρo = get(input, "ρo", 49.1)
    ρw = get(input, "ρw", 64.79)
    #! format: off
    nc, nconn = grid.nc, grid.connlist.nconn
    fluid = OWFluid(nc, nconn, ρo, ρw, pvto, pvtw, krow)
    #! format: on
    return fluid, input
end

function set_init_solution(input::Dict, fluid::AbstractFluid)
    if input["fluid"] == "OW"
        # Update State
        po, sw = input["po"], input["sw"]
        set_fluid_tn(fluid, po, sw)
        update_primary_variable(fluid, po, sw)
    elseif input["fluid"] ∈ ["OIL", "WATER", "GAS"]
        p = input["p"]
        set_fluid_tn(fluid, p)
        update_primary_variable(fluid, p)
    end
end

function Sim(input::Dict)::Sim
    # Init Grid
    (grid, input) = init_grid(input)
    # Init Rock
    (rock, input) = init_rock(input)
    # Construct ConnList
    # println("Construct ConnList")
    construct_connlist(grid, rock)
    sort_conn(grid.connlist)
    # Init State
    (fluid, input) = init_fluid(input, grid)
    # Init Reservoir
    reservoir = StandardReservoir(grid, rock, fluid)
    # Construct Wells
    wells = Dict{String, AbstractFacility}()
    well_names = Set{String}()
    for p in input["producers"]
        name = p["name"]
        @assert !(name in well_names)
        push!(well_names, name)
        wells[name] = init_well(PRODUCER, p, fluid.nv, grid, rock)
    end

    for p in input["injectors"]
        name = p["name"]
        @assert !(name in well_names)
        push!(well_names, name)
        wells[p["name"]] = init_well(INJECTOR, p, fluid.nv, grid, rock)
    end

    # Construct Scheduler
    sch = Scheduler()
    sch.t_current = 0.0
    sch.t0 = 0.0
    sch.dt_max = get(input, "dt_max", 100.)
    sch.dt0 = get(input, "dt0", 0.01)
    reset_dt(sch)
    t_end = get(input, "t_end", 100.)
    set_time_step(sch, get(input, "time_step", [t_end]))

    # Construct Nonlinear Solver
    nsolver = NRSolver()
    nsolver.max_newton_iter = get(input, "max_newton_iter", 10)
    nsolver.min_err = get(input, "min_err", 1.0e-6)
    # Construct Linear Solver
    solver_type = get(input, "linear_solver", "GMRES_ILU0")
    τ = get(input, "τ", 0.1)
    if solver_type == "GMRES_ILU0"
        lsolver = GMRES_ILU0_Solver(τ=τ)
    elseif solver_type == "GMRES_CPR"
        lsolver = GMRES_CPR_Solver(grid.nc, grid.neighbors, τ=τ)
    else
        lsolver = Julia_BackSlash_Solver()
    end

    set_init_solution(input, fluid)
    # println("Update Phases")
    update_phases(fluid, grid.connlist)
    # println("Compute Residual")
    compute_residual(fluid, grid, rock, wells, sch.dt)

    # Init residual and jacobian
    # println("Init residual and jacobian")
    init_nsolver(nsolver, grid, fluid)
    return Sim(reservoir, wells, sch, nsolver, lsolver)
end


function init_well(T::WellType, well_option::Dict, nv::Int, grid::CartGrid, rock::AbstractRock)
    name = well_option["name"]
    perf = Vector{Int}()
    for indices in well_option["perforation"]
        push!(perf, get_grid_index(grid, indices...))
    end
    radius = get(well_option, "radius", 0.5)
    #
    well = StandardWell{T}(name, perf, radius, nv)
    # Set control mode and target
    well.mode = get_ctrl_mode[well_option["mode"]]
    well.target = well_option["target"]
    # Set limits
    limits = get(well_option, "limits", Dict{Limit, Float64}())
    for (k,v) in limits
        well.limits[k] = v
    end
    # Comput well index
    compute_wi(well, grid, rock)
    return well
end

function add_well(sim::Sim, welltype::String, well_option::Dict)::Nothing
    name = well_option["name"]
    @assert !(name in keys(sim.facility))
    T = get_well_type[lowercase(welltype)]
    reservoir = sim.reservoir
    nv, grid, rock = reservoir.fluid.nv, reservoir.grid, reservoir.rock
    sim.facility[name] = init_well(T, well_option, nv, grid, rock)
    reset_dt(sim.scheduler)
    return nothing
end

"""
    change_well_mode(sim::Sim, name::String, mode::String, target::Float64)

Change the control mode of well `name` to be `mode` with target value `target`

# Arguments
- `sim::Sim`: Sim object
- `name::String`: name of the well
- `mode::String`: target mode
    - "bhp": constant BHP
    - "shut": shut-in
    - "orat": constant oil rate
    - "wrat": constant water rate
    - "lrat": constatnt liquid rate
- `target::Float64`

# Examples
```jldoctest
julia> using ResSimAD: get_model, change_well_mode

julia> sim, options = get_model("example1");

julia> change_well_mode(sim, "P1", "bhp", 6200.);

julia> println(sim.facility["P1"].mode)
CBHP

julia> println(sim.facility["P1"].target)
6200.0
```

"""
function change_well_mode(sim::Sim, name::String, mode::String, target::Float64)::Nothing
    @assert name in keys(sim.facility)
    sim.facility[name].mode = get_ctrl_mode[mode]
    sim.facility[name].target = target
    reset_dt(sim.scheduler)
    return nothing
end

"""
    change_well_target(sim::Sim, name::String, target::Float64; cut_dt=false)

Change the control target of well `name` to be `target`. cut time step to be
`sim.scheduler.dt0` if `cut_dt=true`

# Arguments
- `sim::Sim`: Sim object
- `name::String`: name of the well
- `target::Float64`
- `cut_dt::Boolean`

# Examples
```jldoctest
julia> using ResSimAD: get_model, change_well_target

julia> sim, options = get_model("example1");

julia> change_well_target(sim, "P1", 6200.);

julia> println(sim.facility["P1"].target)
6200.0

```

"""
function change_well_target(sim::Sim, name::String, target::Float64; cut_dt=false)::Nothing
    @assert name in keys(sim.facility)
    sim.facility[name].target = target
    if cut_dt
        reset_dt(sim.scheduler)
    end
    return nothing
end

"""
    get_well_rates(sim::Sim, name::String, data::String)

Get column `data` from dataframe `sim.facility[name].results`

# Arguments
- `sim::Sim`: Sim object
- `name::String`: name of the well
- `data::String`: column name
    - "Time": Time steps
    - "ORAT": Oil rate
    - "WRAT": Water rate
    - "GRAT": Gas rate
    - "LRAT": Liquid rate
    - "WBHP": BHP

# Examples
```jldoctest
julia> using ResSimAD: get_model, runsim, get_well_rates, SILENT

julia> sim, options = get_model("example1");

julia> runsim(sim, verbose=SILENT);

julia> tstep = get_well_rates(sim, "P1", "Time");

julia> qo = get_well_rates(sim, "P1", "ORAT");

julia> length(tstep) > 0
true
julia> length(qo) > 0
true
```

"""
function get_well_rates(sim::Sim, name::String, data::String)
    @assert name in keys(sim.facility)
    return sim.facility[name].results[!, data]
end

function shut_well(sim::Sim, name::String)::Nothing
    @assert name in keys(sim.facility)
    sim.facility[name].mode = get_ctrl_mode["shut"]
    return nothing
end


## Define some function for convenience

# Define function po(sim), pw(sim), so(sim), sw(sim) ....
for v in (:p, :s, :b, :μ, :kr, :λ, :ρ, :γ, :ΔΨ, :f, :pn, :sn, :bn, :ρs, :pvt)
    for p in ("o", "w", "g")
        @eval function $(Symbol(v,p))(x::Sim)
            getfield(x, :reservoir).fluid.phases[Symbol($p)].$v
        end
    end
end

# Define function ao(sim), aw(sim), ro(sim), rw(sim) ....
for v in (:a, :r)
    for p in ("o", "w", "g")
        @eval function $(Symbol(v,p))(x::Sim)
            getfield(x, :reservoir).fluid.components[Symbol($p)].$v
        end
    end
end

# Define function po_rec(sim), pw_rec(sim) ...
for v in ("p", "s")
    for p in ("o", "w", "g")
        @eval function $(Symbol(v, p, "_rec"))(x::Sim)
            getfield(getfield(x, :reservoir).fluid.phases[Symbol($p)], Symbol($v, "_rec"))
        end
    end
end

# Define function kx(sim), ϕ(sim) ...
for v in (:kx, :ky, :kz, :ϕ)
    @eval function $v(x::Sim)
        getfield(x, :reservoir).rock.$v
    end
end

# Define function nc(sim), nx(sim) ...
for v in (:nc, :nx, :ny, :nz, :dx, :dy, :dz, :v, :d, :connlist)
    @eval function $v(x::Sim)
        getfield(x, :reservoir).grid.$v
    end
end

# Define function residual, jacobian
for v in (:residual, :jac)
    @eval function $v(x::Sim)
        getfield(x, :nsolver).$v
    end
end

# Define function reservoir(sim), facility(sim) ...
for v in ("reservoir", "facility", "scheduler", "nsolver", "lsolver")
    @eval $(Symbol(v))(x::Sim) = getfield(x, Symbol($v))
end
# Overwrite getproperty function for sim
# sim.$symbol = symbol(sim)
# Now we can use sim.po to get sim.reservoir.fluid.phases.o.p
# Similarly sim.ro => sim.reservoir.fluid.components.o.r
# sim.po_rec => sim.reservoir.fluid.phases.o.p_rec
# At the same time, sim.reservoir sim.facility sim.nsolver ... still works

import Base
Base.getproperty(x::Sim, s::Symbol) = eval(s)(x)
