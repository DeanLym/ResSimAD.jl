using Distributed: addprocs, @everywhere, @spawnat, fetch, workers, rmprocs

# Launch 5 Julia worker processes
nrun = 5
addprocs(nrun);


@everywhere using ResSimAD: get_model, Sim, runsim, get_well_rates

@everywhere function forecast(perm)
    _, options = get_model("example1")
    options["perm"] = perm
    sim = Sim(options);
    runsim(sim);
    t = get_well_rates(sim, "P1", "TIME")
    qo = get_well_rates(sim, "P1", "ORAT")
    qw = get_well_rates(sim, "P1", "WRAT")
    return t, qo, qw
end

perms = rand(nrun) * 200.0 .+ 100.0;

data_ref = []
for (i, worker) in enumerate(workers())
    push!(data_ref, @spawnat worker forecast(perms[i]))
end

results = []
for i = 1:nrun
    push!(results, fetch(data_ref[i]))
end

using Plots
using Plots.PlotMeasures

cmap = cgrad(:jet);

p1 = plot(xlabel="Day", ylabel="Oil rate (STB/Day)", title="P1 oil rate");
p2 = plot(xlabel="Day", ylabel="Water rate (STB/Day)", title="P1 water rate");
for i = 1:nrun
    t, qo, qw = results[i]
    plot!(p1, t, qo, label="run $i", linewidth=2, marker=true)
    plot!(p2, t, qw, label="run $i", linewidth=2, marker=true)
end
plot(p1, p2, layout=(1,2), size=(720, 280), bottom_margin = 10px)

for worker in workers()
    rmprocs(worker)
end
