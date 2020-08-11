# Basic workflow

The basic workflow for ResSimAD.jl consists of three steps:

1. Specify simulation setups
2. Run simulation
3. Visualize simulation results

## Specify simulation setups
The most convenient way for specifying simulation setups for ResSimAD.jl is using
a nested Dictionary. The following example shows how simulation setups are specified
in the [`ResSimAD.Models.example1`](@ref) model.

```@example workflow
using ResSimAD: get_example_data

## Specify input
# Grid and Rock
options = Dict();
options["nx"] = 30; options["ny"] = 15; options["nz"] = 1;
options["dx"] = 50.; options["dy"] = 50.; options["dz"] = 20.;
options["d"] = 8000.;
options["perm"] = 200.; options["poro"] = 0.2;
# Fluid
options["fluid"] = "OW"
options["po"] = 6000.;
options["sw"] = 0.1;
options["PVDO"] = get_example_data("PVDO.DAT");
options["PVTW"] = get_example_data("PVTW.DAT");
options["SWOF"] = get_example_data("SWOF.DAT");
# Wells
options["producers"] = [];
p1 = Dict();
p1["name"] = "P1"; p1["perforation"] = [(1,1,1)]; p1["radius"] = 0.5;
p1["mode"] = "bhp"; p1["target"] = 5500.;
push!(options["producers"], p1);

options["injectors"] = [];
i1 = Dict();
i1["name"] = "I1"; i1["perforation"] = [(30,15,1)]; i1["radius"] = 0.5;
i1["mode"] = "bhp"; i1["target"] = 6500.;
push!(options["injectors"], i1);
# Schedule
options["dt0"] = 0.1
options["dt_max"] = 50.; options["t_end"] = 10 * 182.5;
options["min_err"] = 1.0e-3;

```

Currently supported input options are summarized in [Input options](@ref).

A `Sim` object can then be created with `options`:

```@example workflow
using ResSimAD: Sim

sim = Sim(options);
println("Total number of cells:", sim.nc)
```

## Run simluation
To run simulation, simply call:
```@example workflow
using ResSimAD: runsim

runsim(sim)
```

## Visualize and save simulation results

We use the `Plots` package for results visualizations.

#### Plot production rates
Production rates for each well are stored in a `DataFrame` at
`sim.facility[name].results`, where `name` is the name of the well.
We can also use the [`ResSimAD.get_well_rates`](@ref)
API function to extract production rates.


```@example workflow
using Plots
using Plots.PlotMeasures
using ResSimAD: get_well_rates

gr(format=:svg) # hide
ENV["GKSwstype"] = "100" #hide

# Get time, P1 oil rate, I1 water injection rate
# Option 1
t = sim.facility["P1"].results[!, "TIME"];
qo = sim.facility["P1"].results[!, "ORAT"];
qw = sim.facility["I1"].results[!, "WRAT"];
# Note that column name is case sensitive (all letters must be upper case)

# Option 2
t = get_well_rates(sim, "P1", "TIME");
qo = get_well_rates(sim, "P1", "orat");
qw = get_well_rates(sim, "I1", "Wrat");
# Column name will be converted to upper case inside get_well_rates
# We recommend using Option 2

p1 = plot(t, qo, color=:black, marker=true,
     xlabel="Day", ylabel="Oil rate (STB/Day)",
     title="P1 oil rate");
p2 = plot(t, -qw, color=:black, marker=true,
     xlabel="Day", ylabel="Water inj. Rate (STB/Day)",
     title="I1 water injection Rate");
plot(p1, p2, layout=(1,2), legend=false, size=(720,280), bottom_margin = 10px)
```

#### Plot state maps
Snapshots of the primary variables (`po` and `sw` for this oil-water system) at
each time step are stored in `sim.reservoir.fluid.phases.o.p_rec`
and `sim.reservoir.fluid.phases.w.s_rec`. We can also use derived properties
`sim.po_rec` and `sim.sw_rec` or the API function [`ResSimAD.get_state_map`](@ref)
to access them.

The type for `sim.po_rec` and `sim.sw_rec` is `Dict{Float64, Vector{Float64}}`,
the keys are time step values (`Float64` rounded to `6` digits).

```@example workflow
# Visualize state maps
using ResSimAD: get_grid_index, isproducer, get_state_map
# Get po, sw at the end of simulation
# Option 1
po = sim.reservoir.fluid.phases.o.p_rec[round(t[end], digits=6)];
sw = sim.reservoir.fluid.phases.w.s_rec[round(t[end], digits=6)];
# Option 2
po = sim.po_rec[round(t[end], digits=6)];
sw = sim.sw_rec[round(t[end], digits=6)];
# Option 3
po = get_state_map(sim, "po", t[end]);
sw = get_state_map(sim, "sw", t[end]);
# `t` will be rounded to 6 digits inside get_state_map
# We recommend using Option 3

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
cmap = cgrad(:jet);
p1 = heatmap(reshape(po, sim.nx, sim.ny), color=cmap, title="Po at Day $(t[end])");
plot_wells_2d(p1, sim);
p2 = heatmap(reshape(sw, sim.nx, sim.ny), color=cmap, title="Sw at Day $(t[end])");
plot_wells_2d(p2, sim);
plot(p1, p2, layout=(1,2), size=(700,450))
```

#### Plot newton iterations
Number of newton iterations at each time step is stored in `sim.nsolver.num_iter`,
we can plot this to check the convergence behavior for this simulation run.
```@example workflow
newton_iter = sim.nsolver.num_iter;
scatter(newton_iter, markershape=:square, size=(360,280), legend=false,
        xlabel="Time step", ylabel="Newton iterations")
```
