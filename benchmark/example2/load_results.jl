using HDF5
using ResSimAD
using CSV
using DataFrames
using DelimitedFiles
using Statistics

## Load simulation results
dir = joinpath(pkgdir(ResSimAD), "benchmark", "example2")
results = Dict()
runtimes = Dict()
wnames = ["I1", "I2", "P1", "P2"]

# load ResSimAD results
results["ResSimAD"] = Dict()

data_types = Dict([
    ("P1", [("oil rate", "ORAT"), ("water rate", "WRAT")]),
    ("P2", [("oil rate", "ORAT"), ("water rate", "WRAT")]),
    ("I1", [("water inj. rate", "WRAT")]),
    ("I2", [("water inj. rate", "WRAT")]),
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

indices = Dict([
    ("P1", 1), ("P2", 2), ("I1", 1), ("I2", 2),
])

data_types = Dict([
    ("P1", [("oil rate", "WOPR"), ("water rate", "WWPR")]),
    ("P2", [("oil rate", "WOPR"), ("water rate", "WWPR")]),
    ("I1", [("water inj. rate", "WWIR")]),
    ("I2", [("water inj. rate", "WWIR")])
])

fn = joinpath(dir, "Eclipse", "results", "EXAMPLE2.h5")
fid = h5open(fn, "r")
results["Eclipse"]["Day"] = read(fid, "summary_vectors/TIME/0/values")


for (iw, wn) in enumerate(wnames)
    ii = indices[wn]
    for (key, col) in data_types[wn]
        ind = names(fid["summary_vectors/$col"])[ii]
        results["Eclipse"][wn * " " * key] = read(fid, "summary_vectors/$col/$ind/values")
    end
end
close(fid)

runtimes["Eclipse"] = readdlm(joinpath(dir, "Eclipse", "results", "average_runtime.txt"))

# load MRST results
results["MRST"] = Dict()

data_types = Dict([
    ("P1", [("oil rate", "qo.txt", 3), ("water rate", "qw.txt", 3)]),
    ("P2", [("oil rate", "qo.txt", 4), ("water rate", "qw.txt", 4)]),
    ("I1", [("water inj. rate", "qw.txt", 1)]),
    ("I2", [("water inj. rate", "qw.txt", 2)]),
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
    ("P2", [("oil rate", "OPR"), ("water rate", "WPR",)]),
    ("I1", [("water inj. rate", "WIR")]),
    ("I2", [("water inj. rate", "WIR")])
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
    ("P2", [("oil rate", "WOPR"), ("water rate", "WWPR",)]),
    ("I1", [("water inj. rate", "WWIR")]),
    ("I2", [("water inj. rate", "WWIR")])
])

for (iw, wn) in enumerate(wnames)
    for (key, col) in data_types[wn]
        results["OPM"][wn * " " * key] = read(fid, string(col, ":", wn))
    end
end

close(fid)

runtimes["OPM"] = readdlm(joinpath(dir, "OPM", "results", "average_runtime.txt"))
