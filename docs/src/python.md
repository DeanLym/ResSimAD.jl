# Python usage

With the [PyJulia](https://pyjulia.readthedocs.io/en/latest/) python package, `ResSimAD.jl` is also accessible from `Python`.

## PyJulia setup
Install the `PyJulia` package in python and setup the link between `Julia` and `Python`, following the steps described in [PyJulia](https://pyjulia.readthedocs.io/en/latest/). And make sure that `ResSimAD.jl` is runnable in the `Julia` version that is linked to `Python`.

## Basic workflow in Python
After installing `PyJulia`, calling `import julia` in Python will launch a Julia process in the background. With `from julia import Main`, we can access variable `x` in the background Julia process at `Main.x`.

We can then import the `ResSimAD` module with `from julia import ResSimAD`.

```python
import julia
from julia import Main
from julia import ResSimAD
```

Simulation setups can be specified using a `Python` dictionary.

```python
options = {}
# Grid and Rock
options["nx"] = 30; options["ny"] = 15; options["nz"] = 1;
options["dx"] = 50.; options["dy"] = 50.; options["dz"] = 20.;
options["d"] = 3000.
options["perm"] = 200.; options["poro"] = 0.2;

# Fluid
options["fluid"] = "OW"
options["po"] = 6000.; options["sw"] = 0.1;
options["PVDO"] = ResSimAD.get_example_data("PVDO.DAT");
options["PVTW"] = ResSimAD.get_example_data("PVTW.DAT");
options["SWOF"] = ResSimAD.get_example_data("SWOF.DAT");

# Wells
options["producers"] = [];
p1 = {};
p1["name"] = "P1"; p1["perforation"] = [(1,1,1)]; p1["radius"] = 0.5;
p1["mode"] = "bhp"; p1["target"] = 5500.;
options["producers"].append(p1);

options["injectors"] = [];
i1 = {};
i1["name"] = "I1"; i1["perforation"] = [(30,15,1)]; i1["radius"] = 0.5;
i1["mode"] = "bhp"; i1["target"] = 7000.;
options["injectors"].append(i1);

# Schedule
options["dt0"] = 0.1;
options["dt_max"] = 50.;
options["t_end"] = 1825.0;
options["min_err"] = 1.0e-3
```

Then we can pass the `options` to the `ResSimAD.Sim` function to create a `Sim` object in the background Julia process. This `Sim` object is accessible with `Main.sim`.

```python
Main.sim = ResSimAD.Sim(options)
```

To run simulation, simply call

```python
ResSimAD.runsim(Main.sim)
```

To extract simulation results, run the following code. The returned simulation will be converted to `numpy` arrays. We can then analyze/visualize the results in `Python` with `matplotlib` (not shown here).

```python
t = ResSimAD.get_well_rates(Main.sim, "P1", "TIME")
qo = ResSimAD.get_well_rates(Main.sim, "P1", "ORAT")
qw = ResSimAD.get_well_rates(Main.sim, "P1", "WRAT")
po = ResSimAD.get_state_map(Main.sim, "po", t[-1])
sw = ResSimAD.get_state_map(Main.sim, "sw", t[-1])
```

## Advanced workflow in Python

### Dynamic simulation control in Python
Dynamic simulation control is also doable in `Python`.

```python
import numpy as np
# Create the example1 model
Main.sim, Main.options = ResSimAD.get_model("example1")
# Change well target every 100 days
tsteps = [100., 200., 300.];
for t in tsteps:
    p1_bhp = np.random.rand()*400. + 5500.;
    i1_bhp = np.random.rand()*400. + 6100.;
    ResSimAD.change_well_target(Main.sim, "P1", p1_bhp)
    ResSimAD.change_well_target(Main.sim, "I1", i1_bhp)
    ResSimAD.change_dt(Main.sim, 5.0)
    ResSimAD.step_to(Main.sim, t)
    print(Main.sim.scheduler.t_current)

# Change well target every 5 timesteps
for i in range(2):
    p1_bhp = np.random.rand()*400. + 5500.;
    i1_bhp = np.random.rand()*400. + 6100.;
    ResSimAD.change_well_target(Main.sim, "P1", p1_bhp)
    ResSimAD.change_well_target(Main.sim, "I1", i1_bhp)
    ResSimAD.change_dt(Main.sim, 5.0)
    for j in range(5):
        ResSimAD.time_step(Main.sim)
    print(Main.sim.scheduler.t_current)

# Add new well
i2 = {
    "name": "I2",
    "perforation": [(8,12,1)],
    "radius": 0.5,
    "mode": "bhp",
    "target": 6500.
}

ResSimAD.add_well(Main.sim, "injector", i2)

t_end = 1200.
ResSimAD.step_to(Main.sim, t_end)
print(Main.sim.scheduler.t_current)

```

### Run multiple simulations in parallel in Python

The example in [Run multiple simulations in parallel](@ref) shows how to run multiple simulations in parallel in Julia. The following example shows how to run multiple simulations in Python.

```python
import julia
from julia import Main
```

We can use the `Distributed` Julia module in Python to launch multiple Julia worker processes
```python
from julia import Distributed
Main.nrun = 5
Distributed.addprocs(Main.nrun)
```

Next we need use the `@everywhere`, `@spawnat` macro in Julia. This is a Julia specific syntax, which can not be directly called in Python. But we can use the `Main.eval` method to execute any Julia code in the background Julia process:

```python
Main.eval(
"""
using Distributed

@everywhere using ResSimAD

@everywhere function forecast(perm)
    _, options = get_model("example1")
    options["perm"] = perm
    sim = Sim(options);
    runsim(sim, verbose=SILENT);
    t = get_well_rates(sim, "P1", "TIME")
    qo = get_well_rates(sim, "P1", "ORAT")
    qw = get_well_rates(sim, "P1", "WRAT")
    return t, qo, qw
end

"""
)
```

We can still create variables in Python and pass to Julia:

```python
import numpy as np
Main.perms = np.random.rand(Main.nrun) * 100 + 200
```

Next we launch multiple simulations:

```python
Main.eval("""
data_ref = []
for (i, worker) in enumerate(workers())
    push!(data_ref, @spawnat worker forecast(perms[i]))
end
""")

```

Notice that the `perms` variable in the above Julia code was created in Python and passed to Julia.

Then we run the following code to fetch simulation results from the Julia worker processes to the Julia main process.

```python
Main.eval("""
results = []
for i = 1:nrun
    push!(results, fetch(data_ref[i]))
end
""")
```

The results from the multiple simulations are then accessible in the `Main.results` variable in Python. In this case, the `Main.results` is a `list` of 5 `tuples`, each `tuple` contains three `numpy arrays` - `(t, qo, qw)`

## Use precompiled ResSimAD.jl in Python

We can also use precompiled `ResSimAD.jl` in Python.

The first appraoch is to precompile `ResSimAD.jl` and replace the default Julia system image, as described in [Precompile](@ref).

Alternatively, in PyJulia, we can use a customize system image for Julia, as described in the [Custom Julia system image](https://pyjulia.readthedocs.io/en/latest/sysimage.html) section of the PyJulia documentation. If we have precompiled ResSimAD.jl and created a custom system image called `ResSimADSysImg.so`. Then we call uses this system image to initialize Julia in Python:

```python
from julia import Julia
jl = Julia(sysimage="ResSimADSysImg.so")

from julia import Main
from julia import ResSimAD
```
