# Example models

There are several example models in the `src/models/` folder.
We can use the [`ResSimAD.get_model`](@ref) function to create these example models.

For instance, to create the `example1` model:
```julia
using ResSimAD:get_model, runsim
sim, options = get_model("example1");
runsim(sim);
```

Example code for simulating the example models and plotting results are available in the `examples/` folder.

## List of example models

```@docs
ResSimAD.Models.example1
```

```@docs
ResSimAD.Models.example2
```

```@docs
ResSimAD.Models.example3
```
