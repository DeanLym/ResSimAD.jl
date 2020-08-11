# Advanced workflow

## Dynamic simulation control
`ResSimAD.jl` allows flexible control over simulation runs.

API functions [`ResSimAD.time_step`](@ref) and [`ResSimAD.step_to`](@ref) enable simulating the model for a single timestep, or to a specific time, respectively.

API functions [`ResSimAD.change_well_mode`](@ref), [`ResSimAD.change_well_target`](@ref), [`ResSimAD.shut_well`](@ref) and [`ResSimAD.add_well`](@ref) enable dynamic modification of well setups during a simulation run.

After changing well settings, time step size will be cut down to `sim.scheduler.dt0`. We can also use [`ResSimAD.change_dt`](@ref) to change the time step size.

Here is an example where well control settings are changed dynamically during simulation runs. This is very convenient in the context of well control optimizations.


```@example workflow2
# Create simulation model
using ResSimAD: get_model, change_well_target, change_dt, step_to, time_step

sim, options = get_model("example1");

# Change BHP control every 100 days
tsteps = [100., 200., 300.];
for t in tsteps
    p1_bhp = rand()*400. + 5500.;
    i1_bhp = rand()*400. + 6100.;
    change_well_target(sim, "P1", p1_bhp)
    change_well_target(sim, "I1", i1_bhp)
    change_dt(sim, 5.0)
    step_to(sim, t)
end

# Change BHP every 5 time steps
for i = 1:2
    p1_bhp = rand()*400. + 5500.;
    i1_bhp = rand()*400. + 6100.;
    change_well_target(sim, "P1", p1_bhp)
    change_well_target(sim, "I1", i1_bhp)
    change_dt(sim, 5.0)
    for j = 1:5
        time_step(sim)
    end
end

```

```@example workflow2
## Plot results
using ResSimAD: get_well_rates, get_state_map, get_grid_index, isproducer
using Plots
using Plots.PlotMeasures
cmap = cgrad(:jet);
gr(format=:svg) # hide
ENV["GKSwstype"] = "100" #hide

# BHP
t = get_well_rates(sim, "P1", "TIME");

pw_p1 = get_well_rates(sim, "P1", "WBHP");
pw_i1 = get_well_rates(sim, "I1", "WBHP");
p1 = plot(t, pw_p1, color=:black, marker=true,
     xlabel="Day", ylabel="BHP (psi)",
     title="P1 BHP");
p2 = plot(t, pw_i1, color=:black, marker=true,
     xlabel="Day", ylabel="BHP (psi)",
     title="I1 BHP");
plot(p1, p2, layout=(1,2), legend=false, size=(720,280), bottom_margin = 10px)
```

```@example workflow2
# Oil rate and water injection rate
qo = get_well_rates(sim, "P1", "ORAT");
qw = get_well_rates(sim, "I1", "WRAT");
p1 = plot(t, qo, color=:black, marker=true,
     xlabel="Day", ylabel="Oil rate (STB/Day)",
     title="P1 oil rate");
p2 = plot(t, -qw, color=:black, marker=true,
     xlabel="Day", ylabel="Water inj. rate (STB/Day)",
     title="I1 water injection rate");
plot(p1, p2, layout=(1,2), legend=false, size=(720,280), bottom_margin = 10px)
```

```@example workflow2
# Define a function for plotting wells
function plot_wells_2d(plt, sim; markersize=5, color=:white)
    for w in values(sim.facility)
        i, j, _ = get_grid_index(sim.reservoir.grid, w.ind[1])
        if isproducer(w)
            marker = :circle
        else
            marker = :dtriangle
        end
        scatter!(plt, [j,], [i,], m=(marker, markersize, color), legend=false)
    end
end

# Plot state maps
po = get_state_map(sim, "po", t[end]);
sw = get_state_map(sim, "sw", t[end]);
p1 = heatmap(reshape(po, sim.nx, sim.ny), color=cmap, title="Po at day $(t[end])");
plot_wells_2d(p1, sim);
p2 = heatmap(reshape(sw, sim.nx, sim.ny), color=cmap, title="Sw at day $(t[end])");
plot_wells_2d(p2, sim);
plot(p1, p2, layout=(1,2), size=(700,450))
```

We can also add new well dynamically. This can be useful for field development optimization.

```@example workflow2
using ResSimAD: add_well

i2 = Dict("name" => "I2", "perforation"=>[(8,12,1)],
          "radius"=>0.5, "mode"=>"bhp", "target"=>6500.);

add_well(sim, "injector", i2);

t_end = 1200.;
step_to(sim, t_end);

```

```@example workflow2
# Plot production curves
qo = get_well_rates(sim, "P1", "ORAT");
t2 = get_well_rates(sim, "I2", "TIME");
qw1 = get_well_rates(sim, "I1", "WRAT");
qw2 = get_well_rates(sim, "I2", "WRAT");
p1 = plot(t, qo, color=:black, marker=true, label = "P1",
     xlabel="Day", ylabel="Oil rate (STB/Day)",
     title="P1 oil rate");
p2 = plot(t, -qw1, color=:blue, marker=true, label="I1");
plot!(p2, t2, -qw2, color=:red, marker=true, label="I2",
     xlabel="Day", ylabel="Water inj. rate (STB/Day)",
     title="Water injection rate");
plot(p1, p2, layout=(1,2), size=(720,280), bottom_margin = 10px)
```

```@example workflow2
# Plot state maps
po = get_state_map(sim, "po", t_end);
sw = get_state_map(sim, "sw", t_end);
p1 = heatmap(reshape(po, sim.nx, sim.ny), color=cmap, title="Po at day $(t[end])");
plot_wells_2d(p1, sim);
p2 = heatmap(reshape(sw, sim.nx, sim.ny), color=cmap, title="Sw at day $(t[end])");
plot_wells_2d(p2, sim);
plot(p1, p2, layout=(1,2), size=(700,450))
```


## Run multiple simulations in parallel

With the Julia [Distributed](https://docs.julialang.org/en/v1/manual/distributed-computing/) module, it is very convenient to run multiple simulations in parallel. Here we
provide an example:

First, launch Julia worker processes.

```julia
using Distributed: addprocs, @everywhere, @spawnat, fetch, workers

# Launch 5 Julia worker processes
nrun = 5
addprocs(nrun);

```

Then import `ResSimAD.jl` and define a `forecast` function on all worker processes with the `@everywhere` macro.

```julia
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
    push!(data_ref, @spawnat worker forecast(perms[i]))
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

```@example workflow2
nrun = 5; # hide
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

p1 = plot(xlabel="Day", ylabel="Oil rate (STB/Day)", title="P1 oil rate");
p2 = plot(xlabel="Day", ylabel="Water rate (STB/Day)", title="P1 water rate");
for i = 1:nrun
    t, qo, qw = results[i]
    plot!(p1, t, qo, label="run $i", linewidth=2, marker=true)
    plot!(p2, t, qw, label="run $i", linewidth=2, marker=true)
end
plot(p1, p2, layout=(1,2), size=(720, 280), bottom_margin = 10px)
```
