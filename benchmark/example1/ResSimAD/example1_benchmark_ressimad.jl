using ResSimAD

include(joinpath(pkgdir(ResSimAD), "benchmark", "benchmark_ressimad.jl"))

model_name = "example1"
root_dir = @__DIR__;

run_benchmark_ressimad(model_name, root_dir)
