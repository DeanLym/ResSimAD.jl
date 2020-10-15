
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
    T = get_well_type[lowercase(welltype)]
    reservoir = sim.reservoir
    nv, grid, rock = reservoir.fluid.nv, reservoir.grid, reservoir.rock
    name = well_option["name"]
    if "wi" ∈ keys(well_option)
        parse_well_option(well_option, welltype, well_names, true)
        sim.facility[name] = init_well(T, well_option, nv, grid)
    else # "radius"
        parse_well_option(well_option, welltype, well_names, false)
        sim.facility[name] = init_well(T, well_option, nv, grid, rock)
    end
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


"""
    get_data(sim::Sim, dataname::String)

Get data for dataname `dataname`.

# Examples
```jldoctest
julia> using ResSimAD: get_model, step_to, get_data, silence

julia> silence();

julia> sim, options = get_model("example1");

julia> step_to(sim, 10.0);

julia> λo = get_data(sim, "λo");

julia> size(λo)
(450,)

```
"""
function get_data(sim::Sim, dataname::String)
    ret = @eval $(Symbol(lowercase(dataname)))($sim)
    if typeof(ret) == Array{Float64, 1}
        return ret
    else
        return value(ret)
    end
end

"""
    save_results(sim::Sim; dir::String="./", state_report::String="report")

Save simulation results to `dir`. Well results will be saved to `wellname.csv`.
State maps for primary variables will be saved to `state.h5`.

# Arguments
- `sim::Sim`: Sim object
- `dir::String`: target folder
- `state_report::String`:
    - "nothing": no state maps will be saved
    - "report": state maps at days `sim.scheduler.time_step` will be saved
    - "all": state maps at all time steps will be saved

# Examples
```jldoctest
julia> using ResSimAD: get_model, runsim, silence, save_results

julia> silence();

julia> sim, options = get_model("example1");

julia> runsim(sim);

julia> save_results(sim);

```

"""
function save_results(sim::Sim; dir::String="./", state_report::String="report")
    if !(ispath(dir))
        mkdir(dir)
    end
    info(LOGGER, "Saving simulation results")
    # Save production curves
    for w in values(sim.facility)
        CSV.write(joinpath(dir, w.name * ".csv"), w.results)
    end
    # Save statemaps
    state_report = lowercase(state_report)
    if state_report == "nothing"
        return nothing
    end
    vars = primary_variables(sim.reservoir.fluid)
    if state_report == "report"
        tstep = sim.scheduler.time_step
    elseif state_report == "all"
        tstep = get_well_rates(sim, collect(keys(sim.facility))[1], "time")
    end
    info(LOGGER, "Saving state maps at Days $tstep")
    nt = length(tstep)
    fid = h5open(joinpath(dir, "state.h5"), "w")
    write(fid, "t", tstep)
    for v in vars
        data = zeros(nt, sim.nc)
        for it = 1:nt
            data[it, :] = get_state_map(sim, v, tstep[it])
        end
        write(fid, v, data)
    end
    close(fid)
    return nothing
end
