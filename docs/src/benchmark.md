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
using ResSimAD
using Plots
using Plots.PlotMeasures
gr(format=:svg) # hide
ENV["GKSwstype"] = "100" #hide
# Load results
include(joinpath(pkgdir(ResSimAD), "benchmark", "example1", "load_results.jl"))

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
ResSimAD.jl | GMRES + ILU | 70 | 186
Eclipse | ORTHOMIN + NF  | 77| 178
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
