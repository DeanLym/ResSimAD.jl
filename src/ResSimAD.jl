module ResSimAD

using Memento

include("utils/global.jl")
include("autodiff/autodiff.jl")
include("reservoir/rock/rock.jl")
include("reservoir/grid/grid.jl")
include("reservoir/fluid/fluid.jl")
include("reservoir/reservoir.jl")
include("facility/facility.jl")
include("schedule/schedule.jl")
include("linearsolver/linearsolver.jl")
include("inputparser/inputparser.jl")
include("simmaster/simmaster.jl")
include("models/models.jl")

using .AutoDiff:value

using .Grid: get_grid_index

using .Facility: isproducer

using .SimMaster:Sim, runsim, time_step, step_to, newton_step, add_well, change_well_mode,
            change_well_target, shut_well, get_well_rates, change_dt, get_residual_error,
            get_state_map, set_perm, set_poro


using .Models: get_model, get_example_data

export value
export get_grid_index
export change_dt
export set_perm, set_poro
export Sim, runsim, time_step, step_to, newton_step
export get_model, get_example_data, get_state_map, get_well_rates
export get_residual_error
export change_well_mode, change_well_target, shut_well, add_well

function __init__()
    logger = getlogger()

    logger.handlers["console"] = DefaultHandler(
        stdout, DefaultFormatter("[{level}]: {msg}"),
        Dict{Symbol, Any}(
            :colors => Dict{AbstractString, Symbol}(
                "trace" => :normal,
                "debug" => :blue,
                "info" => :normal,
                "notice" => :cyan,
                "warn" => :magenta,
                "error" => :red,
                "critical" => :yellow,
                "alert" => :white,
                "emergency" => :red,
            )
        )
    )
end

"""
    silence()

Do not show any message except for "error" messages. See `ResSimAD.verbose()` for more information.

# Examples
```jldoctest
julia> using ResSimAD:silence

julia> silence();

```

"""
function silence()
    setlevel!(getlogger("root"), "error"; recursive=true)
end

"""
    debug()

Show "debug" messages. See `ResSimAD.verbose()` for more information.

# Examples
```jldoctest
julia> using ResSimAD:debug

julia> debug();

```

"""
function debug()
    setlevel!(getlogger("root"), "debug"; recursive=true)
end

"""
    verbose(level::String)

Set verbose level. Show messages with higher level than the set `level`.
See Memento.jl for more information.

# Arguments
- `level::String`: verbose levels (low to high)
    - `"debug"`: verbose message used for debugging
    - `"info"`: general information
    - `"notice"`: important events that are part of normal execution
    - `"warn"`: warning that may cause the program to fail
    - `"error"`: error that causes the program to terminate

# Examples
```jldoctest
julia> using ResSimAD:verbose

julia> verbose("debug");

```

"""
function verbose(level::String)
    setlevel!(getlogger("root"), level; recursive=true);
end

end # module
