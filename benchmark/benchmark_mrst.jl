using ResSimAD
using Memento
using Statistics
using Dates
using DelimitedFiles

function run_benchmark_mrst(model_name, dir)
    cd(dir)
    ## Setup logging
    # Get system information
    machine = Sys.MACHINE
    # Add time stamp
    time_now = now()
    timestamp = Dates.format(now(), "dd_u_yyyy_HH_MM_SS")
    logger = getlogger("root")
    log_file = joinpath(dir, "mrst_$(model_name)_$(machine)_$(timestamp).log")
    push!(logger, DefaultHandler(log_file,  DefaultFormatter("[{level}]: {msg}")))

    logs = []
    runtimes = []
    for irun = 1:5
        info(logger, "MRST run $irun")

        # Launch MRST
        out = Pipe();
        cmd1 = `matlab.exe -nodesktop -nosplash -nojvm -minimize -r "run ('example1.m'), quit" -wait -log`
        run(pipeline(cmd1, stderr=out));
        close(out.in);

        # Read runtime from log
        log = String(read(out));
        secs = match(r"\d+", match(r"control steps.*Seconds", log).match).match
        millisecs = match(r"\d+", match(r"\d* Milliseconds \*\*\*", log).match).match
        elapsed_time = parse(Float64, secs) + parse(Float64, millisecs) / 1000.
        push!(runtimes, elapsed_time)
        push!(logs, log)
        info(logger, "Elapsed time for run $(irun): $(round(runtimes[irun], digits=3)) seconds\n")
    end

    info(logger, "Average run time for the 5 simulations: $(round(mean(runtimes), digits=3)) seconds")

    # Write run logs to log file
    for (irun, log) in enumerate(logs)
        info(logger, "\nRun $irun log:")
        info(logger, log)
    end
    # Remove log-file from handlers
    logger.handlers = filter!(x -> x[1]=="console", logger.handlers)

    writedlm(joinpath(dir, "average_runtime.txt"), runtimes)

    results_dir = joinpath(dir, "results")
    if !(isdir(results_dir))
        mkdir(results_dir)
    end

    files = ("qo.txt", "qw.txt", "t.txt", "average_runtime.txt")
    for file in files
        mv(joinpath(dir,file), joinpath(results_dir, file); force=true)
    end


end
