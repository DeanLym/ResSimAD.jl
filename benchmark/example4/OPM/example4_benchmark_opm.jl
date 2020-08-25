using ResSimAD
"""
Simulate example4 with OPM Flow
"""

include(joinpath(pkgdir(ResSimAD), "benchmark", "benchmark_opm.jl"))

model_name = "example4"

root_dir = @__DIR__;

run_benchmark_opm(model_name, root_dir; nrun=3)


