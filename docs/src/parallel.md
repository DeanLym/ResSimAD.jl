# Multiple simulations in parallel

With the Julia [Distributed](https://docs.julialang.org/en/v1/manual/distributed-computing/) module, it is very convenient to run multiple simulations in parallel. Here we
provide an example:

First, launch Julia worker processes.
```@example parallel
using Distributed: addprocs, @everywhere, @spawnat, fetch, workers

# Launch 5 Julia worker processes
nrun = 5
addprocs(nrun)

# Check number of workers
println("Number of workers:", length(workers()))
```

Then import `ResSimAD.jl` and define a `forecase` function on all worker processes with the `@everywhere` macro.

```@example parallel
@everywhere using ResSimAD: get_model, Sim, runsim, SILENT, get_well_rates

@everywhere function forecast(perm)
    _, options = get_model("example1")
    options["perm"] = perm
    sim = Sim(options);
    runsim(sim, verbose=SILENT);
    t = get_well_rates(sim, "P1", "TIME")
    qo = get_well_rates(sim, "P1", "ORAT")
    qw = get_well_rates(sim, "P1", "WRAT")
    return t, qo, qw
end

```

Next, launch simulations on the worker processes with the `@spawnat` macro.

```julia
perms = rand(nrun) * 200.0 .+ 100.0;

data_ref = []
for (i, worker) in enumerate(workers())
    push!(data_ref, @spawnat worker forecast(perm[i]))
end

```

Next, retrieve simulation results with the `fetch` function

```julia
results = []
for i = 1:nrun
    push!(results, fetch(data_ref[i]))
end
```

We can then plot simulation results

```@example parallel
perms = rand(nrun) * 200.0 .+ 100.0; # hide
using ResSimAD # hide
function forecast(perm) # hide
    _, options = get_model("example1") # hide
    options["perm"] = perm # hide
    sim = Sim(options); # hide
    runsim(sim, verbose=SILENT); # hide
    t = get_well_rates(sim, "P1", "TIME") # hide
    qo = get_well_rates(sim, "P1", "ORAT") # hide
    qw = get_well_rates(sim, "P1", "WRAT") # hide
    return t, qo, qw # hide
end # hide
results = [] # hide
for i = 1:nrun push!(results, forecast(perms[i])) end # hide
using Plots
using Plots.PlotMeasures

gr(format=:svg) # hide
ENV["GKSwstype"] = "100" #hide

p1 = plot(xlabel="Day", ylabel="Oil rate (STB/Day)", title="P1 oil rate");
p2 = plot(xlabel="Day", ylabel="Water rate (STB/Day)", title="P1 water rate");
for i = 1:nrun
    t, qo, qw = results[i]
    plot!(p1, t, qo, label="run $i", linewidth=2, marker=true)
    plot!(p2, t, qw, label="run $i", linewidth=2, marker=true)
end
plot(p1, p2, layout=(1,2), size=(720, 280), bottom_margin = 10px)
```
