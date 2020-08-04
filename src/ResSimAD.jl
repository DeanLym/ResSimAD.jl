module ResSimAD

include("utils/global.jl")
include("utils/inputparse.jl")
include("autodiff/autodiff.jl")
include("reservoir/rock/rock.jl")
include("reservoir/grid/grid.jl")
include("reservoir/fluid/fluid.jl")
include("reservoir/reservoir.jl")
include("facility/facility.jl")
include("schedule/schedule.jl")
include("linearsolver/linearsolver.jl")
include("simmaster/simmaster.jl")

using .AutoDiff:value

using .SimMaster:Sim, runsim, step, step_to, newton_step, add_well, change_well_mode,
            change_well_target, shut_well, get_well_rates

export value
export Sim, runsim, step, step_to, newton_step

end # module
