println("Importing ResSimAD (1st time) takes:")
@time using ResSimAD

println("Creating Sim object (1st time) takes:")
@time sim, options = get_model("example1")

println("Running simulation (1st time) takes:")

@time runsim(sim; verbose=SILENT);


println("\nImporting ResSimAD (2nd time) takes:")
@time using ResSimAD

println("Creating Sim object (2nd time) takes:")
@time sim, options = get_model("example1")

println("Running simulation (2nd time) takes:")

@time runsim(sim; verbose=SILENT);


using ResSimAD

include(joinpath(pkgdir(ResSimAD), "examples", "precompile.jl"))
