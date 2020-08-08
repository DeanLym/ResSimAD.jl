using Revise
using ResSimAD: get_model, get_well_rates, get_state_map
using ResSimAD: change_well_target, add_well, step_to, time_step, change_dt
using Plots
using Plots.PlotMeasures
cmap = cgrad(:jet);

include("utils.jl");
## Create simulation model

sim, options = get_model("example1");

# Change BHP control every 100 days
tsteps = [100., 200., 300.];

for t in tsteps
    p1_bhp = rand()*400. + 5500.;
    i1_bhp = rand()*400. + 6100.;
    change_well_target(sim, "P1", p1_bhp);
    change_well_target(sim, "I1", i1_bhp);
    change_dt(sim, 5.0);
    step_to(sim, t);
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


## Plot results

# production curves
t = get_well_rates(sim, "P1", "TIME");

pw_p1 = get_well_rates(sim, "P1", "WBHP");
pw_i1 = get_well_rates(sim, "I1", "WBHP");
p1 = plot(t, pw_p1, color=:black, marker=true,
     xlabel="Day", ylabel="BHP (psi)",
     title="P1 BHP");
p2 = plot(t, pw_i1, color=:black, marker=true,
     xlabel="Day", ylabel="BHP (psi)",
     title="I1 BHP");
plot(p1, p2, layout=(1,2), legend=false, size=(720,300), bottom_margin = 10px)

qo = get_well_rates(sim, "P1", "ORAT");
qw = get_well_rates(sim, "I1", "WRAT");
p1 = plot(t, qo, color=:black, marker=true,
     xlabel="Day", ylabel="Oil rate (STB/Day)",
     title="P1 oil rate");
p2 = plot(t, -qw, color=:black, marker=true,
     xlabel="Day", ylabel="Water inj. rate (STB/Day)",
     title="I1 water injection rate");
plot(p1, p2, layout=(1,2), legend=false, size=(720,280), bottom_margin = 10px)


# Plot state maps
po = get_state_map(sim, "po", t[end]);
sw = get_state_map(sim, "sw", t[end]);
p1 = heatmap(reshape(po, sim.nx, sim.ny), color=cmap, title="Po at day $(t[end])");
plot_wells_2d(p1, sim);
p2 = heatmap(reshape(sw, sim.nx, sim.ny), color=cmap, title="Sw at day $(t[end])");
plot_wells_2d(p2, sim);
plot(p1, p2, layout=(1,2), size=(700,450))


## Add new wells

i2 = Dict("name" => "I2", "perforation"=>[(8,12,1)],
          "radius"=>0.5, "mode"=>"bhp", "target"=>6500.);

add_well(sim, "injector", i2);

t_end = t[end] + 500.
step_to(sim, t_end);

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
     title="water injection rate");
plot(p1, p2, layout=(1,2), size=(720,280), bottom_margin = 10px)

# Plot state maps
po = get_state_map(sim, "po", t_end);
sw = get_state_map(sim, "sw", t_end);
p1 = heatmap(reshape(po, sim.nx, sim.ny), color=cmap, title="Po at day $(t[end])");
plot_wells_2d(p1, sim);
p2 = heatmap(reshape(sw, sim.nx, sim.ny), color=cmap, title="Sw at day $(t[end])");
plot_wells_2d(p2, sim);
plot(p1, p2, layout=(1,2), size=(700,450))
