module ResSimAD

include("AutoDiff.jl")
include("Grid.jl")
include("State.jl")
include("Schedule.jl")
include("Solver.jl")
include("SimCtrl.jl")


using .AutoDiff:param, zeros_tensor, Tensor
using .State:OWState
using .SimCtrl:Sim, setup, runsim

export param, zeros_tensor, Tensor
export Sim, setup, runsim

# function setnv(nv::Int)
#     global Nv = nv
# end

end # module
