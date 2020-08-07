# Input Options

```@meta
DocTestSetup = quote
    options = Dict()
    Nc = 2
end
```

---

## Grid
---
### dx, dy, dz
Cell size along $x$, $y$, $z$ dimension
- Data type

Type    | Size
:---:   | :---:
Float64 | -
Vector{Float64} | $N_c$

- Examples:

```jldoctest
julia> options["dx"] = 50.;

julia> options["dy"] = 50.;

julia> options["dz"] = 10.;

```

```jldoctest
julia> options["dx"] = 50. * ones(Nc);

julia> options["dy"] = 50. * ones(Nc);

julia> options["dz"] = 10. * ones(Nc);

```
---
### d
Cell center depths
- Data type

Type    | Size
:---:   | :---:
Float64 | -
Vector{Float64} | $N_c$

- Example:

```jldoctest
julia> options["d"] = 8000.;

```

```jldoctest
julia> options["d"] = 8000. * ones(Nc);

```
---
### nx, ny, nz
Number of cells along $x$, $y$, $z$ dimension

- Data type

Type    | Size
:---:   | :---:
Int64 | -

- Example:

```jldoctest
julia> options["nx"] = 50;

julia> options["ny"] = 60;

julia> options["nz"] = 1;

```
---
## Rock
---
### perm
Cell permeability value

- Data type

Type    | Size
:---:   | :---:
Float64 | -
Vector{Float64} | $N_c$

- Example:

```jldoctest
julia> options["perm"] = 200.;

```

```jldoctest
julia> options["perm"] = 200. * ones(Nc);

```
---
### poro
Cell porosity value

- Data type

Type    | Size
:---:   | :---:
Float64 | -
Vector{Float64} | $N_c$

- Example:

```jldoctest
julia> options["poro"] = 0.2;

```

```jldoctest
julia> options["poro"] = 0.2 * ones(Nc);

```
---
## Fluid
---
### fluid
Fluid system type.

- Data type

Type    | Size
:---:   | :---:
String | -

- Supported values
    - "OW" - oil-water system

- Examples

```jldoctest
julia> options["fluid"] = "OW";

```
---
### po
Initial oil phase pressure

- Data type

Type    | Size
:---:   | :---:
Float64 | -
Vector{Float64} | $N_c$

- Examples

```jldoctest
julia> options["po"] = 6000.0;

```

```jldoctest
julia> options["po"] = 6000.0 * ones(Nc);

```
---
### sw
Initial water phase saturation

- Data type

Type    | Size
:---:   | :---:
Float64 | -
Vector{Float64} | $N_c$

- Examples

```jldoctest
julia> options["sw"] = 0.1;

```

```jldoctest
julia> options["sw"] = 0.1 * ones(Nc);

```
---
### PVDO
File name for the dead oil PVT table

- Data type

Type    | Size
:---:   | :---:
String | -

```jldoctest
julia> options["PVDO"] = "PVDO.DAT";

```
---
### PVTW
File name for the water PVT function

- Data type

Type    | Size
:---:   | :---:
String | -

```jldoctest
julia> options["PVTW"] = "PVTW.DAT";

```
---
### SWOF
File name for the oil-water two-phase relative permeabilty table

- Data type

Type    | Size
:---:   | :---:
String | -

```jldoctest
julia> options["SWOF"] = "SWOF.DAT";

```

---
## Wells
---
### producers, injectors
List of producers/injectors

- Data type

Type    | Size
:---:   | :---:
Vector{Dict} | number of producers/injectors

- Each dictionary contains the key-value pairs
    - `"name" -> String` - well name
    - `"perforation" -> Vector{Tuple{Int, Int, Int}}` - indices for perforated grid blocks
    - `"radius" -> Float64` - wellbore radius
    - `"mode" -> String` - well control mode
        - `"bhp"`: constant BHP
        - `"shut"`: shut-in
        - `"orat"`: constant oil rate
        - `"wrat"`: constant water rate
        - `"lrat"`: constatnt liquid rate
    - `"target" -> Float64` - well control target

- Example

```jldoctest
julia> options["producers"] = [];

julia> p1 = Dict();

julia> p1["name"] = "P1"; p1["perforation"] = [(1,1,1)]; p1["radius"] = 0.5;

julia> p1["mode"] = "bhp"; p1["target"] = 5500.;

julia> push!(options["producers"], p1);

```

```jldoctest
julia> options["injectors"] = [];

julia> i1 = Dict();

julia> i1["name"] = "I1"; i1["perforation"] = [(30,15,1)]; i1["radius"] = 0.5;

julia> i1["mode"] = "bhp"; i1["target"] = 6500.;

julia> push!(options["injectors"], i1);

```

---
## Schedule
---
### dt0
Initial time step size

- Data type

Type    | Size
:---:   | :---:
Float64 | -

- Example

```jldoctest
julia> options["dt0"] = 0.1;

```
---
### dt_max
Maximum time step size

- Data type

Type    | Size
:---:   | :---:
Float64 | -

- Example

```jldoctest
julia> options["dt_max"] = 50.0;

```

---
### t_end
Maximum time step size

- Data type

Type    | Size
:---:   | :---:
Float64 | -

- Example

```jldoctest
julia> options["t_end"] = 1825.0;

```
---
## Nonlinear solver
---
### min_err
Minimum residual error for convergence

- Data type

Type    | Size
:---:   | :---:
Float64 | -

- Example

```jldoctest
julia> options["min_err"] = 1.0e-3;

```
