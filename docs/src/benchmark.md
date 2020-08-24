# Benchmark

In this section, `ResSimAD.jl` is compared with `Eclipse`, `ADGPRS`, `MRST` and `OPM`.

Input files, scripts, results and output logs for all simulators are available in `ResSimAD.jl/benchmark/`. 
To reproduce the benchmarking results, proper installation, environmental variables setup and valid licenses for `Eclipse`, `ADGPRS`, `MRST` and `OPM` are required.

The version, machine and CPU information for producing the benchmarking results presented here is summarized in the table below.

Simulator    | Version | Machine | CPU
:---:   | :---: | :---: | :---:
ResSimAD.jl | dev | Windows 10 (x86_64-w64-mingw32)| Intel(R) Core(TM) i7-4790 CPU @ 3.60GHz
Eclipse | v2017.2 | Windows 10 (x86_64-w64-mingw32)| Intel(R) Core(TM) i7-4790 CPU @ 3.60GHz
MRST | 2020a | Windows 10 (x86_64-w64-mingw32)| Intel(R) Core(TM) i7-4790 CPU @ 3.60GHz
ADGPRS | 1.9 | Ubuntu-20.04 on WSL2 (x86_64-pc-linux-gnu)| Intel(R) Core(TM) i9-9920X CPU @ 3.50GHz
OPM | 2020.04 | Ubuntu-20.04 on WSL2 (x86_64-pc-linux-gnu)| Intel(R) Core(TM) i9-9920X CPU @ 3.50GHz

Note that at the time these benchmarking results were generated, we ended up using two different machines, due to practical reasons such as different supported platforms from 
different simulators, availability of licenses (`Eclipse` and matlab license for `MRST`). Therefore, the performances of the simulators
are subject to small differences. These differences are in general small and will not affect the conclusions presented in this benchmarking study.

We thank the Stanford Energy Resources Engineering department and the Stanford Smart Fields Consortium for providing `ADGPRS`, license for `Eclipse` and license for Matlab.

## Example1

The `example1` model is a 2D deal-oil model (30x15 grid) with one injector and one producer. 
This simple model is a good candidate for testing the overhead, especially those introduced by AD (except for `Eclipse`, all other simulators use AD).

### Simluation results comparison

Simulation results for `example1` from all simulators are stored in `ResSimAD.jl/benchmark/example1`. The following code reads and plots these stored simulation results.

```@example benchmark
using HDF5
using ResSimAD
using Plots
using Plots.PlotMeasures
using CSV
using DataFrames
using DelimitedFiles
using Statistics

gr(format=:svg) # hide
ENV["GKSwstype"] = "100" #hide

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


```



```@example benchmark
# Plot production rates
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
```

### Performance comparison

Each simulator was ran for 5 times. The average run times are plotted below. 
Run times for `ResSimAD.jl` were directly measured. 
Run times for other simulators were parsed from the output logs. 
This way, the overhead of launching the simulators was not included, especially for `MRST` where launching Matlab took quite some time.
For all simulators, the default linear solver was used. File IO for `Eclipse`, `OPM` and `ADGPRS` was minimized as much as possible.


Simulator    | Linear Solver | Time steps | Newton Iterations
:---:   | :---: | :---: | :---:
ResSimAD.jl | GMRES + ILU | 70 | 168
Eclipse | ORTHOMIN + NF  | 77| 197
MRST | Matlab Backslash | 70 | 176
ADGPRS | GMRES + CPR0 | 72 | 197
OPM | BiCG-stab + ILU0 | 68 | 224

```@example benchmark
# Plot average run time
p1 = plot(ylabel="Average run time (seconds)", legend=false);
for key in sort!(collect(keys(runtimes)))
    bar!(p1, [key], [mean(runtimes[key])], color=:gray)
end

p2 = plot(ylabel="Average run time (seconds)", legend=false);
for key in sort!(["Eclipse", "ResSimAD", "ADGPRS", "OPM"])
    bar!(p2, [key], [mean(runtimes[key])], color=:gray)
end

plot(p1, p2, layout=(1,2), size=(720, 280), bottom_margin = 10px)

```
