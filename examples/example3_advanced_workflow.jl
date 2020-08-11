using Revise
using ResSimAD

using Plots
using Plots.PlotMeasures

include("utils.jl");

sim, options = get_model("example3");

runsim(sim)

# production curves
producers = ["P1", "P2", "P3"]
columns = ["WBHP", "ORAT", "WRAT", "LRAT", ]
ylabels = ["BHP (psi)", "Oil rate (STB/Day)", "Water rate (STB/Day)", "Liquid rate (STB/Day)"]
titles = ["BHP", "Oil rate", "Water rate", "Liquid rate"]
plt = []
for (i, col) in enumerate(columns)
    p = plot(xlabel="Day", ylabel=ylabels[i], title=titles[i])
    for (j, w) in enumerate(producers)
        x = get_well_rates(sim, w, "TIME");
        y = get_well_rates(sim, w, col);
        plot!(p, x, y, label=w, marker=true)
    end
    push!(plt, p)
end
plot(plt..., layout=(2,2), size=(720,600), legend = true)

injectors = ["I1", "I2"]
columns = ["WBHP", "WRAT"]
ylabels = ["BHP (psi)", "Water inj. rate (STB/Day)"]
titles = ["BHP", "Water injection rate"]
plt = []
for (i, col) in enumerate(columns)
    p = plot(xlabel="Day", ylabel=ylabels[i], title = titles[i])
    for (j, w) in enumerate(injectors)
        x = get_well_rates(sim, w, "TIME");
        y = get_well_rates(sim, w, col);
        plot!(p, x, abs.(y), label=w, marker=true)
    end
    push!(plt, p)
end

plot(plt..., layout=(1,2), size=(720,300), legend = true)



# Plot state maps
cmap = cgrad(:jet);
t = sim.scheduler.t_current
po = get_state_map(sim, "po", t);
sw = get_state_map(sim, "sw", t);
p1 = heatmap(reshape(po, sim.nx, sim.ny), color=cmap, title="Po at day $(t)");
plot_wells_2d(p1, sim);
p2 = heatmap(reshape(sw, sim.nx, sim.ny), color=cmap, title="Sw at day $(t)");
plot_wells_2d(p2, sim);
plot(p1, p2, layout=(1,2), size=(700,280))
