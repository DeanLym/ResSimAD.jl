using HDF5
using ResSimAD
using Plots
using Plots.PlotMeasures
using CSV
using DataFrames
using DelimitedFiles
using Statistics

## Load simulation results
dir = joinpath(pkgdir(ResSimAD), "benchmark", "example1")
results = Dict()
runtimes = Dict()
wnames = ["I1", "P1"]

# load ResSimAD results
results["ResSimAD"] = Dict()

data_types = Dict([
    ("P1", [("oil rate", "ORAT"), ("water rate", "WRAT")]),
    ("I1", [("water inj. rate", "WRAT")])
])

for (iw, wn) in enumerate(wnames)
    fn = joinpath(dir, "ResSimAD", "results", "$wn.csv")
    df = DataFrame(CSV.File(fn))
    results["ResSimAD"]["Day"] = df[!, "TIME"]
    for (key, col) in data_types[wn]
        results["ResSimAD"][wn * " " * key] = df[!, col]
    end
end

runtimes["ResSimAD"] = readdlm(joinpath(dir, "ResSimAD", "results", "average_runtime.txt"))
# Load Eclipse results

results["Eclipse"] = Dict()

data_types = Dict([
    ("P1", [("oil rate", "WOPR"), ("water rate", "WWPR")]),
    ("I1", [("water inj. rate", "WWIR")])
])

fn = joinpath(dir, "Eclipse", "results", "EXAMPLE1.h5")
fid = h5open(fn, "r")
results["Eclipse"]["Day"] = read(fid, "summary_vectors/TIME/0/values")

for (iw, wn) in enumerate(wnames)
    for (key, col) in data_types[wn]
        ind = names(fid["summary_vectors/$col"])[1]
        results["Eclipse"][wn * " " * key] = read(fid, "summary_vectors/$col/$ind/values")
    end
end
close(fid)

runtimes["Eclipse"] = readdlm(joinpath(dir, "Eclipse", "results", "average_runtime.txt"))

# load MRST results
results["MRST"] = Dict()

data_types = Dict([
    ("P1", [("oil rate", "qo.txt", 2), ("water rate", "qw.txt", 2)]),
    ("I1", [("water inj. rate", "qw.txt", 1)])
])

results["MRST"]["Day"] = readdlm(joinpath(dir, "MRST", "results", "t.txt"), ',')

for (iw, wn) in enumerate(wnames)
    for (key, file, ind) in data_types[wn]
        data = readdlm(joinpath(dir, "MRST", "results", file), ',')[:, ind]
        results["MRST"][wn * " " * key] = data
    end
end

runtimes["MRST"] = readdlm(joinpath(dir, "MRST", "results", "average_runtime.txt"))

# load ADGPRS resuls
results["ADGPRS"] = Dict()

file = joinpath(dir, "ADGPRS", "results", "OUTPUT.rates.txt")
header = String.(readdlm(file)[1,:])
data = readdlm(file,skipstart=1)
df = DataFrame(data)
rename!(df, header)

data_types = Dict([
    ("P1", [("oil rate", "OPR"), ("water rate", "WPR",)]),
    ("I1", [("water inj. rate", "WIR")])
])

results["ADGPRS"]["Day"] = df[!, "Day"]

for (iw, wn) in enumerate(wnames)
    for (key, col) in data_types[wn]
        results["ADGPRS"][wn * " " * key] = df[!, string(wn, ":", col)] *6.28981
    end
end

runtimes["ADGPRS"] = readdlm(joinpath(dir, "ADGPRS", "results", "average_runtime.txt"))

# load OPM results
results["OPM"] = Dict()

file = joinpath(dir, "OPM", "results", "results.h5")
fid = h5open(file, "r")

results["OPM"]["Day"] = read(fid, "YEARS") * 365.

data_types = Dict([
    ("P1", [("oil rate", "WOPR"), ("water rate", "WWPR",)]),
    ("I1", [("water inj. rate", "WWIR")])
])

for (iw, wn) in enumerate(wnames)
    for (key, col) in data_types[wn]
        results["OPM"][wn * " " * key] = read(fid, string(col, ":", wn))
    end
end

close(fid)

runtimes["OPM"] = readdlm(joinpath(dir, "OPM", "results", "average_runtime.txt"))


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
to_plot = ["P1 oil rate", "P1 water rate", "I1 water inj. rate"]
for key in to_plot
    p = plot(xlabel="Day", ylabel=key * " (stb/day)", size=(420, 320),
             legend=:right, title=key)
    for case in sort(collect(keys(results)))
        plot!(p, results[case]["Day"][7:3:end], abs.(results[case][key][7:3:end]), label=case,
                line=(linetypes[case], 3.0), marker=markers[case])
    end
    push!(plts, p)
end

plot(plts..., layout=(3,1), size=(320,820), left_margin=20px)


p1 = plot(ylabel="Average run time (seconds)", legend=false);
for key in sort!(collect(keys(runtimes)))
    bar!(p1, [key], [mean(runtimes[key])], color=:gray)
end

p2 = plot(ylabel="Average run time (seconds)", legend=false);
for key in sort!(["Eclipse", "ResSimAD", "ADGPRS", "OPM"])
    bar!(p2, [key], [mean(runtimes[key])], color=:gray)
end

plot(p1, p2, layout=(1,2), size=(720, 280), bottom_margin = 10px)
