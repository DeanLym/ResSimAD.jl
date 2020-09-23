# ResSimAD.jl

[![Build Status](https://travis-ci.com/DeanLym/ResSimAD.jl.svg?token=zPX8pK8q8xHrqbTxACjW&branch=master)](https://travis-ci.com/DeanLym/ResSimAD.jl)

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://deanlym.github.io/ResSimAD.jl/dev/)

A light-weighted reservoir simulator in a modern programming language.

## Features
- **Interactivity**: run simulations in interactive programming environments such as `Jupyter notebook` and `VSCode`. With the `PyJulia` python module, `ResSimAD.jl` is also accessible from `Python`.
- **Dynamic control**: run simulation timestep by timestep, or newton step by newton step, and modify parameters or controls dynamically.
- **Distributed computing**: with the `Distributed` Julia built-in module, it is very convenient to run multiple simulations in parallel. In addition, with the ClusterManagers.jl package, it is very convenient to run multiple simulations in parallel on a HPC cluster. This is very useful for tasks such as field optimization and history matching that require simulating a large number of models.

## Automatic Differentiation (AD)
An operator-overloading-based forward mode AD framework is developed to compute gradients automatically. Instead of using existing Julia AD packages such as [Zygote.jl](https://github.com/FluxML/Zygote.jl), [Tracker.jl](https://github.com/FluxML/Tracker.jl) or [DualNumbers.jl](https://github.com/JuliaDiff/DualNumbers.jl), a customized AD framework is developed here for the sake of efficiency. The AD framework underlying `ResSimAD.jl` is tailored for the operations in reservoir simulation. This allows maximum level of optimization which makes this AD framework almost as fast as hand-written differentiation. On the other hand, the AD framework may need to be extended if some new operations are introduced when extending the functionality of `ResSimAD.jl`.

## Functionality
Functionality-wise, `ResSimAD.jl` is still at a very early stage. It currently works for simple simulation models:
- Grid: 3D Cartesian grid
- Fluid: two phase (oil-water) dead oil
- Well: single perforation

But the underlying framework of `ResSimAD.jl` is designed for easy extension to more complex simulation models, such as those with unstructured grid, three-phase black oil fluid model, wells with multiple perforations. This is largely facilitated by the powerful type system and multiple dispatch in Julia, and the AD framework.

## Performance
For the supported cases, the speed of `ResSimAD.jl` is comparable to professional simulators written in low-level languages including [Eclipse](https://www.software.slb.com/products/eclipse), [OPM](https://opm-project.org/) and [ADGPRS](https://supri-b.stanford.edu/research-areas/ad-gprs). See [`Benchmark`] section in the documentation for benchmark comparisons with [MRST](https://www.sintef.no/projectweb/mrst/), [ADGPRS](https://supri-b.stanford.edu/research-areas/ad-gprs), [Eclipse](https://www.software.slb.com/products/eclipse) and [OPM](https://opm-project.org/).

## Installation
`ResSimAD.jl` is being tested internally. You will need to get access to the repository.

To install, first open `Julia` and press key `]` to activate the package manager.

Then, install the `duneistl_bicgilu_jll.jl` package with

```julia
add https://github.com/DeanLym/duneistl_bicgilu_jll.jl.git
```

Next, install `ResSimAD.jl` with

```julia
add https://github.com/DeanLym/ResSimAD.jl.git
```

The `duneistl_bicgiul_jll.jl` is a Julia wrapper for the BiGCstab solver and the incomplete LU preconditioner from the [DUNE-ISTL](https://dune-project.org/) library. It is the fastest linear solver (2-3x speedup than the default solver) among all the supported linear solvers in `ResSimAD.jl`. 

Set `options["linear_solver"]="BICGSTAB_ILU_DUNE_ISTL"` to use this solver. 

Normally, when installing `ResSimAD.jl`, all dependencies will be automatically installed. However, the `duneistl_bicgiul_jll.jl` has not been registered to the `Julia` registry. So unfortunately we need to install this one manually for now.

## Documentation

Here is the [documentation](https://deanlym.github.io/ResSimAD.jl/dev/).

