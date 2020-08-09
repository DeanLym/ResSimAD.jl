# Precompile


## Why precompiling

Quote from the [PackageCompiler.jl](https://julialang.github.io/PackageCompiler.jl/dev/)
documentation:

> Julia is, in general, a "just-barely-ahead-of-time" compiled language. When you call a function for the first time, Julia compiles it for precisely the types of the arguments given. This can take some time. All subsequent calls within that same session use this fast compiled function, but if you restart Julia you lose all the compiled work.

In a new Julia session, it may feel slow when importing ResSimAD, creating `Sim`
object, running simulation for the first time. This is because there are a lot of
functions getting compiled "just-barely-ahead-of-time". These include functions
in the ResSimAD.jl package itself and many methods in the dependencies of ResSimAD.jl.

The following code demonstrates this:

```julia
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

```

There are huge differences in the execution time between the first calls and the second calls:

```julia
Importing ResSimAD (1st time) takes:
  2.962693 seconds (6.16 M allocations: 384.068 MiB)
Creating Sim object (1st time) takes:
 28.540179 seconds (59.03 M allocations: 2.926 GiB, 5.07% gc time)
Running simulation (1st time) takes:
  1.797259 seconds (4.20 M allocations: 364.556 MiB, 5.76% gc time)

Importing ResSimAD (2nd time) takes:
  1.004753 seconds (2.59 M allocations: 122.898 MiB, 1.96% gc time)
Creating Sim object (2nd time) takes:
  0.001944 seconds (2.95 k allocations: 1.226 MiB)
Running simulation (2nd time) takes:
  0.373982 seconds (138.59 k allocations: 160.657 MiB, 7.20% gc time)
```

Creating and simulating this small model (with only `450 cells`) for the first time will take about `30 secs`. And all subsequent runs will take less than `0.5 secs`. In a long running Julia session (e.g., in Juno or Jupyter notebook), the slow start up time is often times acceptable. However, there are cases where this long start up become unacceptable. For such cases, it is recommended to use the [PackageCompiler.jl](https://julialang.github.io/PackageCompiler.jl/dev/) package to speed up the start up time.

## How to precompile

To precompile ResSimAD.jl, first install the PackageCompiler.jl package. Then run:

```julia
using ResSimAD

fn = joinpath(pkgdir(ResSimAD), "precompile", "precompile.jl");

include(fn)

```

This will take several minutes. You will see the following outputs:

```julia
Precompiling ResSimAD.jl.
Default system image will be replaced.
[ Info: PackageCompiler: creating system image object file, this might take a while...
[ Info: PackageCompiler: default sysimg replaced, restart Julia for the new sysimg to be in effect

```

After the precompiling finishes, restart Julia and run the following code again

```julia
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

```

This time, importing ResSimAD.jl, creating `Sim` object and running simulation for the first time becomes almost as fast as the second calls:

```julia
Importing ResSimAD (1st time) takes:
  0.001425 seconds (1.16 k allocations: 58.562 KiB)
Creating Sim object (1st time) takes:
  0.028607 seconds (63.63 k allocations: 4.854 MiB)
Running simulation (1st time) takes:
  0.399527 seconds (193.98 k allocations: 163.970 MiB, 4.50% gc time)

Importing ResSimAD (2nd time) takes:
  0.000831 seconds (751 allocations: 39.719 KiB)
Creating Sim object (2nd time) takes:
  0.002119 seconds (3.04 k allocations: 1.227 MiB)
Running simulation (2nd time) takes:
  0.352242 seconds (138.59 k allocations: 160.657 MiB, 4.27% gc time)

```




## What precompiling does

To do this, we first create a precompile execution file (available at `precompile/precompule_execution_file.jl`). This file executes the most commonly used ResSimAD.jl API functions. Next, copy this file to some folder and run the following code,
```julia
using PackageCompiler
create_sysimage(:ResSimAD; sysimage_path="ResSimADSysImg.so",
                precompile_execution_file="precompile_execution_file.jl")
```
a custom system image `ResSimADSysImg` will be created in

 When we provide this file to the package compiler, all functions and subroutines that get executed in this file will be compiled. The compiled functions will be stored in a custom system image. We can specify Jul

The same applies to plotting production curves and state maps. The first line plot
and the first heat map take much longer to show up than later plots. This is often
referred to as the "time to first plot" problem.
