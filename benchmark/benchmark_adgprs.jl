using ResSimAD
using Memento
using Statistics
using Dates
using DelimitedFiles

function run_benchmark_adgprs(model_name, dir)
    cd(dir)
    ## Setup logging
    # Get system information
    machine = Sys.MACHINE
    # Add time stamp
    time_now = now()
    timestamp = Dates.format(now(), "dd_u_yyyy_HH_MM_SS")
    logger = getlogger("root")
    log_file = joinpath(dir, "adgprs_$(model_name)_$(machine)_$(timestamp).log")
    push!(logger, DefaultHandler(log_file,  DefaultFormatter("[{level}]: {msg}")))

    # Add timestamp
    info(logger, string(time_now))
    # Add system information to the log_file
    cpu_info = Sys.cpu_info()
    info(logger, "Sys.MACHINE: $(machine)")
    info(logger, "CPU model: $(cpu_info[1].model)")
    info(logger, "Number of CPU cores: $(length(cpu_info))")
    info(logger, "Total memory $(round(Sys.total_memory() / 2^30, digits=1)) GB")
    info(logger, "Free memory $(round(Sys.free_memory() / 2^30, digits=1)) GB\n")

    ## Benchmark ADGPRS run
    logs = []
    runtimes = []
    for irun = 1:5
        info(logger, "ADGPRS run $irun")
        # Launch ADGPRS, direct log to a pipe
        out = Pipe();
        cmd1 = `ADGPRS $(uppercase(model_name)*".DATA") 1 0`;
        run(pipeline(cmd1, stdout=out));
        close(out.in);

        # Read elapsed time from log
        log = String(read(out));
        ret = match(r"\d*\.\d* seconds", log);
        elapsed_time = parse(Float64, match(r"\d*\.\d*", ret.match).match);
        push!(runtimes, elapsed_time)
        push!(logs, log)
        info(logger, "Elapsed time for run $(irun): $(round(runtimes[irun], digits=3)) seconds\n")
    end

    info(logger, "Average run time for the 5 simulations: $(round(mean(runtimes), digits=3)) seconds\n")

    # Write run logs to log file
    for (irun, log) in enumerate(logs)
        info(logger, "\nRun $irun log:")
        info(logger, log)
    end
    # Remove log-file from handlers
    logger.handlers = filter!(x -> x[1]=="console", logger.handlers)

    writedlm(joinpath(dir, "average_runtime.txt"), runtimes)

    ## Remove some files
    files = readdir()

    for file in files
        if endswith(file, ".log.txt") || contains(file, "OUTPUT.FIPS.txt") || file == "fort.66"
            rm(file)
        end
    end

    results_dir = joinpath(dir, "results")
    if !(isdir(results_dir))
        mkdir(results_dir)
    end

    files = ("OUTPUT.rates.txt", "average_runtime.txt")
    for file in files
        mv(joinpath(dir,file), joinpath(results_dir, file); force=true)
    end

end
