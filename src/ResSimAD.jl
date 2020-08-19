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
            get_state_map, set_perm, set_poro, save_results

using .Models: get_model, get_example_data

include("utils/log.jl")

export value
export get_grid_index
export change_dt
export set_perm, set_poro
export Sim, runsim, time_step, step_to, newton_step
export get_model, get_example_data, get_state_map, get_well_rates, save_results
export get_residual_error
export change_well_mode, change_well_target, shut_well, add_well
export add_log_file


end # module
