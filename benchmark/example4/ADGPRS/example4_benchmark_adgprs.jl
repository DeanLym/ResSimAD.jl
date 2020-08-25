using ResSimAD
"""
Simulate example4 with ADGPR
"""

include(joinpath(pkgdir(ResSimAD), "benchmark", "benchmark_adgprs.jl"))

model_name = "example4"
root_dir = @__DIR__;

run_benchmark_adgprs(model_name, root_dir; nrun=3)
