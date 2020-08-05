# Basic workflow

The basic workflow for ResSimAD.jl consists of three steps:

1. Specify simulation setups
2. Run simulation
3. Visualize and save simulation results

## Specify simulation setups
The most convenient way for specifying simulation setups for ResSimAD.jl is using
a nested Dictionary. The following example shows how simulation setups are specified
in the [`ResSimAD.Models.example1`](@ref) model.

```@example workflow
using ResSimAD
# This example will use the following functions/types from ResSimAD
# get_table, runsim, Sim, BRIEF, get_well_rates


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
options["PVDO"] = get_table("PVDO.DAT");
options["PVTW"] = get_table("PVTW.DAT");
options["SWOF"] = get_table("SWOF.DAT");
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

Currently supported input options are summarized in [Input Options](@ref).

A `Sim` object can then be created with `options`:

```@example workflow
sim = Sim(options);
println("Total number of cells:", sim.nc)
```

## Run simluation
To run simulation, simply call:
```@example workflow
runsim(sim; verbose=BRIEF)
```

## Visualize and save simulation results

We use the `Plots` package for results visualizations.

#### Plot production curve
```@example workflow
using Plots
gr(format=:svg) # hide
ENV["GKSwstype"] = "100" #hide

t = get_well_rates(sim, "P1", "TIME");
qo = get_well_rates(sim, "P1", "ORAT");
qw = get_well_rates(sim, "I1", "WRAT");
p1 = plot(t, qo, color=:black, marker=true,
     xlabel="Day", ylabel="Oil Rate (STB/Day)",
     title="P1 Oil Rate");
p2 = plot(t, -qw, color=:black, marker=true,
     xlabel="Day", ylabel="Water Inj. Rate (STB/Day)",
     title="I1 Water Injection Rate");
plot(p1, p2, layout=(1,2), legend=false, size=(720,280))
```

#### Plot state maps
```@example workflow
# Visualize state maps
cmap = cgrad(:jet);
po = reshape(sim.po_rec[t[end]], sim.nx, sim.ny);
sw = reshape(sim.sw_rec[t[end]], sim.nx, sim.ny);
p1 = heatmap(po, color=cmap, title="Po at Day $(t[end])");
p2 = heatmap(sw, color=cmap, title="Sw at Day $(t[end])");
plot(p1, p2, layout=(1,2), size=(720,280))
```

#### Plot newton iterations
```@example workflow
newton_iter = sim.nsolver.num_iter;
scatter(newton_iter, markershape=:square, size=(360,280), legend=false,
        xlabel="Time step", ylabel="Newton iterations")
```
