using ResSimAD
"""
Simulate example1 with OPM Flow
"""

include(joinpath(pkgdir(ResSimAD), "benchmark", "benchmark_opm.jl"))

model_name = "example1"

root_dir = @__DIR__;

log = run_benchmark_opm(model_name, root_dir)
ret = match(r"Total time.*", log)

ret