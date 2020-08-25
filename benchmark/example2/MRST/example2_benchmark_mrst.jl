using ResSimAD
"""
Simulate example2 with MRST
"""

include(joinpath(pkgdir(ResSimAD), "benchmark", "benchmark_mrst.jl"))

model_name = "example2"
root_dir = @__DIR__;

run_benchmark_mrst(model_name, root_dir)
