using ResSimAD
using Memento
using Statistics

model_name = "example2"
"""
The first run includes compilation time.
Simulation time for the first run is not logged.
"""
sim, options = get_model(model_name)
runsim(sim)

"""
Log run time for the second run
"""
# Get system information
machine = Sys.MACHINE
# Add a log file
ResSimAD.add_log_file(joinpath(@__DIR__, "$model_name_$machine.log"))
# Add system information to the log_file
logger = getlogger("root")
cpu_info = Sys.cpu_info()
info(logger, "Sys.MACHINE: $(machine)")
info(logger, "CPU model: $(cpu_info[1].model)")
info(logger, "Number of CPU cores: $(length(cpu_info))")
info(logger, "Total memory $(round(Sys.total_memory() / 2^30, digits=1)) GB")
info(logger, "Free memory $(round(Sys.free_memory() / 2^30, digits=1)) GB\n")

# Run simulation for 5 times and log runtime
runtimes = []
for irun = 1:5
    info(logger, "Simulation run $(irun)")
    t0 = time()

    sim, options = get_model(model_name)
    runsim(sim)

    push!(runtimes, time() - t0)
    info(logger, "Elapsed time for run $(irun): $(round(runtimes[irun], digits=3)) seconds\n")
end

info(logger, "Average run time for the 10 simulations: $(round(mean(runtimes), digits=3)) seconds")
