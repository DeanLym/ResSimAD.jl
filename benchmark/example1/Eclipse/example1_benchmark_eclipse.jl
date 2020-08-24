using ResSimAD
"""
Simulate example1 with Eclipse
"""

include(joinpath(pkgdir(ResSimAD), "benchmark", "benchmark_eclipse.jl"))

model_name = "example1"
root_dir = @__DIR__;

run_benchmark_eclipse(model_name, root_dir)
