using Plots
using Plots.PlotMeasures

ENV["GKSwstype"] = "100" #hide

include(joinpath(@__DIR__, "load_results.jl"))

## Plot results
linetypes = Dict([
    ("ResSimAD", :path),
    ("MRST", :scatter),
    ("Eclipse", :scatter),
    ("ADGPRS", :scatter),
    ("OPM", :scatter),
])

markers = Dict([
    ("ResSimAD", :hline),
    ("MRST", :diamond),
    ("Eclipse", :circle),
    ("ADGPRS", :hexagon),
    ("OPM", :star4),
])

plts = []
to_plot = ["P1 oil rate", "P1 water rate",
            "P2 oil rate", "P2 water rate",
            "I1 water inj. rate", "I2 water inj. rate"]

for key in to_plot
    p = plot(xlabel="Day", ylabel=key * " (stb/day)", size=(420, 320),
            legend=:left, title=key)
    for case in sort(collect(keys(results)))
    # for case in ["Eclipse", "ResSimAD", "ADGPRS", "OPM"]
        plot!(p, results[case]["Day"][3:2:end], abs.(results[case][key][3:2:end]), label=case,
                line=(linetypes[case], 3.0), marker=markers[case])
    end
    push!(plts, p)
end

plot(plts..., layout=(3,2), size=(640,820), left_margin=20px, bottom_margin=10px)


# Plot average run time
p1 = plot(ylabel="Average run time (minutes)", legend=false);
for key in sort!(collect(keys(runtimes)))
    bar!(p1, [key], [mean(runtimes[key]) / 60.], color=:gray)
end

plot(p1, size=(360, 280), bottom_margin = 10px)