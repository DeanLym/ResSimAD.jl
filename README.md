# ResSimAD.jl

[![Build Status](https://travis-ci.com/DeanLym/ResSimAD.jl.svg?token=zPX8pK8q8xHrqbTxACjW&branch=master)](https://travis-ci.com/DeanLym/ResSimAD.jl)

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://deanlym.github.io/ResSimAD.jl/dev/)

A reservoir simulator in a modern programming language.

## Features
- Interactivity: run simulations in interactive programming environments such as Jupyter notebook or Juno. With the PyJulia python module, ResSimAD.jl is also accessible from Python.
- Dynamic control: run simulation timestep by timestep, or newton step by newton step, and modify parameters or controls dynamically.
- Real-time feedback: extract reservoir dynamic states for analysis and visualization.
- Distributed computing: with packages such as ClusterManagers.jl, it is very convenient to run multiple simulations in parallel on a HPC cluster. This is very useful for field optimization and history matching.

## Automatic Differentaition (AD)
An operator-overloading-based forward mode AD framework is developed to compute gradients automatically. Instead of using existing Julia AD packages such as Zygote, Tracker or DualNumbers, a customized AD framework is developed here for the sake of efficiency. The AD framework underlying ResSimAD.jl is tailored for the operations in reservoir simulation. This allows maximum level of optimization which makes this AD framework almost as fast as hand-written differentiation. On the other hand, the AD framework may need to be extended if some new operations are introduced when extending the functionality of ResSimAD.jl.

## Functionality
Functionality-wise, ResSimAD.jl is still at a very early stage. It currently works for simple simulation models:
- Grid: 3D Cartesian grid
- Fluid: two phase (oil-water) dead oil
- Well: single perforation

But the underlying framework of ResSimAD.jl is designed for easy extension for more complex simulation models, such as those with unstructured grid, three-phase black oil fluid model, wells with multiple perforations. This is largely facilitated by the powerful type system and multiple dispatch in Julia, and the AD framework.

Extensions to more complex problems, such as compositional fluid models, multi-segment wells, require more effort. The current AD framework is highly optimized for black oil models and standard wells. Extending the AD framework to handle compositional models, advanced wells and maintaining the same level of efficiency will require some effort. These will not be the focus for ResSimAD.jl for now.

## Performance
For a limited number of test problems, the speed of ResSimAD.jl is comparable to, and in some cases even faster than, Eclipse and ADGPRS. This is very impressive given that ResSimAD.jl is written in pure Julia and that ResSimAD.jl uses AD. This is largely facilitated by
- The fast speed of the Julia language, especially the highly optimized array broadcasting, which makes it possible to completely avoid the overhead of memory allocation induced by operator-overloading-based AD.
- ResSimAD.jl currently uses the GMRES + ILU linear solver from the IterativeSolver.jl package. For the simple problems that ResSimAD.jl currently supports, it is efficient enough.
- File IO can be completely avoid when using ResSimAD.jl in an interactive environment. Simulation results can stay in memory before they are analyzed or visualized for downstream tasks.

## Installation
ResSimAD.jl is currently a private repository. If you have access to this repository,
please download the source code. Then install manually with

```julia
] dev ResSimAD.jl
```

## Examples
Example models are available in src/models/
Code for simulating the example models is available in examples/

## Documentation
Please see this link for [documentation](https://deanlym.github.io/ResSimAD.jl/dev/).
