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
using .SimCtrl:Sim, runsim
using .Well:StandardWell

export param, zeros_tensor, Tensor
export Sim, runsim
export StandardWell

end # module
