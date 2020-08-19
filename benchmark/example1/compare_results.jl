using HDF5
using ResSimAD
using Plots
using CSV
using DataFrames
## Load ResSimAD Results
dir = joinpath(pkgdir(ResSimAD), "benchmark", "example1")

wnames = ["I1", "P1"]
data_types = Dict([
    ("P1", [("ORAT", "WOPR"), ("WRAT", "WWPR")]),
    ("I1", [("WRAT", "WWIR")])
])
titles = Dict([
    ("WOPR", " oil prod. rate (STB/Day)"),
    ("WWPR", " water prod. rate (STB/Day)"),
    ("WWIR", " inj. rate (STB/Day)")
])


fn2 = joinpath(dir, "Eclipse", "results", "EXAMPLE1.h5")
fid = h5open(fn2, "r")
t2 = read(fid, "summary_vectors/TIME/0/values")

plts = []
for (iw, wn) in enumerate(wnames)
    fn1 = joinpath(dir, "ResSimAD", "results", "$wn.csv")

    data1 = DataFrame(CSV.File(fn1))
    for (data1_type, data2_type) in data_types[wn]
        p = plot(xlabel="Day", ylabel=wn * titles[data2_type],
                 size=(420, 280), legend=:right)
        # Plot ResSimAD results
        plot!(p, data1[!, "TIME"], abs.(data1[!, data1_type]), marker=true,
                label="ResSimAD", linewidth=3.0)
        ind = names(fid["summary_vectors/$data2_type"])[1]
        data2 = read(fid, "summary_vectors/$data2_type/$ind/values")
        plot!(p, t2[2:end], data2[2:end], label="Eclipse",
            linewidth=3.0)
        push!(plts, p)
    end
end
close(fid)

for p in plts
    display(p)
end
