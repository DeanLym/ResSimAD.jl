using Revise
using ResSimAD:get_model, runsim, get_well_rates, get_state_map
using Plots

include("utils.jl");

## Create simulation model
sim, options = get_model("example2");

## Run simluation
runsim(sim)

## Plot results
prods, injs = ["P1", "P2"], ["I1", "I2"]
# production curves
plt = []
for w in prods
    t = get_well_rates(sim, w, "TIME");
    qo = get_well_rates(sim, w, "ORAT");
    qw = get_well_rates(sim, w, "WRAT");
    push!(plt, plot(t, qo, color=:black, marker=true,
         xlabel="Day", ylabel="Oil rate (STB/Day)",
         title="$p oil rate"))
    push!(plt, plot(t, qw, color=:black, marker=true,
          xlabel="Day", ylabel="Water rate (STB/Day)",
          title="$w water rate"))
end
plot(plt..., layout=(2,2), legend=false, size=(720,580))


plt = []
for w in injs
    t = get_well_rates(sim, w, "TIME");
    qw = get_well_rates(sim, w, "WRAT");
    push!(plt, plot(t, -qw, color=:black, marker=true,
          xlabel="Day", ylabel="Water inj. rate (STB/Day)",
          title="$w water inj. rate"))
end
plot(plt..., layout=(1,2), legend=false, size=(720,290))


# Plot logk
logk = log.(sim.kx);
p = heatmap(reshape(logk, sim.nx, sim.ny), color=cmap, title="log-permeability");
plot_wells_2d(p, sim);
plot(p, size=(300, 280))

# state maps
cmap = cgrad(:jet);
po = get_state_map(sim, "po", t[end]);
sw = get_state_map(sim, "sw", t[end]);

p1 = heatmap(reshape(po, sim.nx, sim.ny), color=cmap, title="Po at Day $(t[end])");
plot_wells_2d(p1, sim);
p2 = heatmap(reshape(sw, sim.nx, sim.ny), color=cmap, title="Sw at Day $(t[end])");
plot_wells_2d(p2, sim);

plot(p1, p2, layout=(1,2), size=(720,280))

# Newton iterations
newton_iter = sim.nsolver.num_iter;
scatter(newton_iter, markershape=:square, size=(360,280), legend=false,
        xlabel="Time step", ylabel="Newton iterations")
