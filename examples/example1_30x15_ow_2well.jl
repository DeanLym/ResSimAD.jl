using Revise
using ResSimAD
using Plots

sim, options = ResSimAD.get_model("example1")

ResSimAD.runsim(sim)

sim.facility["P1"].results

## Plot results

#
# println("Saving Simulation Results")
# for p in values(sim.producers)
#     CSV.write("$(p.name).txt", p.results)
# end
# for p in values(sim.injectors)
#     CSV.write("$(p.name).txt", p.results)
# end
#
# using DataFrames
#
#
# CSV.write("Pend.txt", DataFrame([sim.state.p_rec[end]], [:Po]))
#
# CSV.write("Swend.txt", DataFrame([sim.state.sw_rec[end]], [:Sw]))
