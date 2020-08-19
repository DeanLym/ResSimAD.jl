using ResSimAD
using Memento
using Statistics
using Dates
using DelimitedFiles

function run_benchmark_ressimad(model_name, dir)
    # Suppress ResSimAD log information
    logger = getlogger("root")
    setlevel!(getlogger("ResSimAD"), "error"; recursive=true)
    """
    The first run includes compilation time.
    Simulation time for the first run is not logged.
    """
    sim, options = get_model(model_name)
    runsim(sim)
    results_dir = joinpath(dir, "results")
    save_results(sim; dir=results_dir)
    """
    Log run time from the second run
    """
    # Get system information
    machine = Sys.MACHINE
    # Add time stamp
    time_now = now()
    timestamp = Dates.format(now(), "dd_u_yyyy_HH_MM_SS")
    # Add a log file
    ResSimAD.add_log_file(joinpath(dir, "ressimad_$(model_name)_$(machine)_$(timestamp).log"))
    # Add timestamp
    info(logger, string(time_now))
    # Add system information to the log_file
    cpu_info = Sys.cpu_info()
    info(logger, "Sys.MACHINE: $(machine)")
    info(logger, "CPU model: $(cpu_info[1].model)")
    info(logger, "Number of CPU cores: $(length(cpu_info))")
    info(logger, "Total memory $(round(Sys.total_memory() / 2^30, digits=1)) GB")
    info(logger, "Free memory $(round(Sys.free_memory() / 2^30, digits=1)) GB\n")

    # Add sim model information
    info(logger, string(sim))

    # Run simulation for 5 times and log runtime
    runtimes = []
    for irun = 1:5
        info(logger, "ResSimAD simulation run $(irun)")
        t0 = time()

        sim, options = get_model(model_name)
        runsim(sim)

        push!(runtimes, time() - t0)
        info(logger, "Elapsed time for run $(irun): $(round(runtimes[irun], digits=3)) seconds\n")
    end

    info(logger, "Average run time for the 5 simulations: $(round(mean(runtimes), digits=3)) seconds")

    # Save average run time to file
    writedlm(joinpath(results_dir, "average_runtime.txt"), runtimes)
end
