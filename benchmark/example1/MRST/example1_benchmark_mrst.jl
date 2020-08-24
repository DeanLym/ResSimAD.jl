using ResSimAD
"""
Simulate example1 with Eclipse (v2017.2)
To run this script, add Eclipse bin folder to the PATH environmental variable
"""

include(joinpath(pkgdir(ResSimAD), "benchmark", "benchmark_mrst.jl"))

model_name = "example1"
root_dir = @__DIR__;

run_benchmark_mrst(model_name, root_dir)
