using Revise
using ResSimAD:get_model, runsim, get_well_rates, get_state_map
using Plots

## Create Simulation model
sim, options = get_model("example1");

## Run Simluation
runsim(sim)

## Plot results

# production curves
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

# state maps
cmap = cgrad(:jet);
po = get_state_map(sim, "po", t[end]);
sw = get_state_map(sim, "sw", t[end]);
p1 = heatmap(reshape(po, sim.nx, sim.ny), color=cmap, title="Po at Day $(t[end])");
p2 = heatmap(reshape(sw, sim.nx, sim.ny), color=cmap, title="Sw at Day $(t[end])");
plot(p1, p2, layout=(1,2), size=(720,280))

# Newton iterations
newton_iter = sim.nsolver.num_iter;
scatter(newton_iter, markershape=:square, size=(360,280), legend=false,
        xlabel="Time step", ylabel="Newton iterations")
