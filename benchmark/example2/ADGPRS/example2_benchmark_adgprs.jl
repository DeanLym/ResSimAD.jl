using ResSimAD
"""
Simulate example1 with ADGPR (v2017.2)
"""

include(joinpath(pkgdir(ResSimAD), "benchmark", "benchmark_adgprs.jl"))

model_name = "example2"
root_dir = @__DIR__;

run_benchmark_adgprs(model_name, root_dir)
