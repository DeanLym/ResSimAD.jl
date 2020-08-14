
function compute_residual(sim::Sim)
    reservoir = sim.reservoir
    grid, fluid, rock = reservoir.grid, reservoir.fluid, reservoir.rock
    facility,sch = sim.facility, sim.scheduler
    compute_residual(fluid, grid, rock, facility, sch.dt)
end


"""
    add_well(sim::Sim, welltype::String, well_option::Dict)

Add one well of type `welltype`, with well setups specified in `well_option`

# Arguments
- `sim::Sim`: Sim object
- `welltype::String`: type of the well (case insensitive)
    - "injector": injector
    - "producer": producer
- `well_option::Dict`: well setups with keys
    - "name" -> String, name of the well
    - "perforation" -> Vector{Tuple{Int, Int, Int}}, indices for perforated blocks
    - "radius" -> Float64, well radius
    - "mode" -> String, well control mode
    - "target" -> Float64, well constrol target

# Examples
```jldoctest
julia> using ResSimAD: get_model, add_well, silence

julia> silence();

julia> sim, options = get_model("example1");

julia> length(keys(sim.facility))
2

julia> p2 = Dict("name" => "P2", "perforation"=>[(5,5,1)], "radius"=>0.5, "mode"=>"bhp", "target"=>5600.);

julia> add_well(sim, "producer", p2);

julia> length(keys(sim.facility))
3

```
"""
function add_well(sim::Sim, welltype::String, well_option::Dict)::Nothing
    well_names = Set{String}(keys(sim.facility))
    parse_well_option(well_option, welltype, well_names)
    T = get_well_type[lowercase(welltype)]
    reservoir = sim.reservoir
    nv, grid, rock = reservoir.fluid.nv, reservoir.grid, reservoir.rock
    name = well_option["name"]
    sim.facility[name] = init_well(T, well_option, nv, grid, rock)
    reset_dt(sim.scheduler)
    sim.nsolver.recompute_residual = true
    notice(LOGGER, "Added new $(welltype) $(name)")
    return nothing
end

"""
    change_well_mode(sim::Sim, name::String, mode::String, target::Float64)

Change the control mode of well `name` to be `mode` with target value `target`

# Arguments
- `sim::Sim`: Sim object
- `name::String`: name of the well
- `mode::String`: target mode, not case sensitive
    - producers:
        - "bhp": constant BHP
        - "shut": shut-in
        - "orat": constant oil rate
        - "wrat": constant water rate
        - "lrat": constatnt liquid rate
    - injectors:
        - "bhp": constant BHP
        - "wrat": constant water rate
- `target::Float64`

# Examples
```jldoctest
julia> using ResSimAD: get_model, change_well_mode, silence

julia> silence();

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
    sim.facility[name].mode = get_ctrl_mode[lowercase(mode)]
    sim.facility[name].target = target
    reset_dt(sim.scheduler)
    sim.nsolver.recompute_residual = true
    notice(LOGGER, "Well $(name) changed to $(sim.facility[name].mode) at $(round(target, digits=3))")
    return nothing
end

"""
    change_well_target(sim::Sim, name::String, target::Float64)

Change the control target of well `name` to be `target`.

# Arguments
- `sim::Sim`: Sim object
- `name::String`: name of the well
- `target::Float64`: well control target

# Examples
```jldoctest
julia> using ResSimAD: get_model, change_well_target, silence

julia> silence();

julia> sim, options = get_model("example1");

julia> change_well_target(sim, "P1", 6200.);

julia> println(sim.facility["P1"].target)
6200.0

```

"""
function change_well_target(sim::Sim, name::String, target::Float64)::Nothing
    @assert name in keys(sim.facility)
    sim.facility[name].target = target
    reset_dt(sim.scheduler)
    sim.nsolver.recompute_residual = true
    notice(LOGGER, "Well $(name) changed to $(sim.facility[name].mode) at $(round(target, digits=3))")
    return nothing
end


"""
    shut_well(sim::Sim, name::String)

Shut well `name`

"""
function shut_well(sim::Sim, name::String)::Nothing
    @assert name in keys(sim.facility)
    sim.facility[name].mode = get_ctrl_mode["shut"]
    reset_dt(sim.scheduler)
    sim.nsolver.recompute_residual = true
    notice(LOGGER, "Well $(name) changed to $(sim.facility[name].mode)")
    return nothing
end

"""
    change_dt(sim::Sim, dt::Float64)

Set next time step size to be `dt`
"""
function change_dt(sim::Sim, dt::Float64)
    set_dt(sim.scheduler, dt)
    sim.nsolver.recompute_residual = true
    notice(LOGGER, "Time step changed to $(round(dt, digits=3))")
end

"""
    get_well_rates(sim::Sim, name::String, data::String)

Get column `data` from dataframe `sim.facility[name].results`

# Arguments
- `sim::Sim`: Sim object
- `name::String`: name of the well
- `data::String`: column name (case insensitive)
    - "TIME": Time steps
    - "ORAT": Oil rate
    - "WRAT": Water rate
    - "GRAT": Gas rate
    - "LRAT": Liquid rate
    - "WBHP": BHP

# Examples
```jldoctest
julia> using ResSimAD: get_model, runsim, get_well_rates, silence

julia> silence();

julia> sim, options = get_model("example1");

julia> runsim(sim);

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
    return sim.facility[name].results[!, uppercase(data)]
end

"""
    get_state_map(sim::Sim, var::String, t::Float64)

Get state map for variable `var` at time `t`.

# Examples
```jldoctest
julia> using ResSimAD: get_model, runsim, get_state_map, get_well_rates, silence

julia> silence();

julia> sim, options = get_model("example1");

julia> runsim(sim);

julia> t = get_well_rates(sim, "P1", "TIME");

julia> po = get_state_map(sim, "po", t[end]);

julia> size(po)
(450,)

```
"""
function get_state_map(sim::Sim, var::String, t::Float64)
    return @eval $(Symbol(lowercase(var)*"_rec"))($sim)[round($t, digits=6)]
end
