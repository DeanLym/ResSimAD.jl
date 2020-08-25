using ResSimAD
"""
Simulate example4 with Eclipse
"""

include(joinpath(pkgdir(ResSimAD), "benchmark", "benchmark_eclipse.jl"))

model_name = "example4"
root_dir = @__DIR__;

run_benchmark_eclipse(model_name, root_dir; nrun=3)
