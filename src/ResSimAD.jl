module ResSimAD

include("Global.jl")
include("AutoDiff.jl")
include("Grid.jl")
include("Fluid.jl")
include("State.jl")
include("Well.jl")
include("Schedule.jl")
include("Solver.jl")
include("SimCtrl.jl")

using .AutoDiff:param, zeros_tensor, Tensor
using .State:OWState
using .Well:StandardWell

using .SimCtrl:Sim, runsim, step, step_to, newton_step, add_well, change_well_mode,
            change_well_target, shut_well, get_well_rates

export param, zeros_tensor, Tensor
export Sim, runsim, step, step_to, newton_step
export StandardWell

end # module
