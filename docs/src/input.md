# Input options

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
String | -

- Example:

```jldoctest
julia> options["perm"] = 200.;

```

```jldoctest
julia> options["perm"] = 200. * ones(Nc);

```

```jldoctest
julia> options["perm"] = "PERM.DAT";

```

Example input file format

```

PERM
-- lines starting with "--" are comments
-- empty lines before keyword are allowed
-- empty lines after / are allowed
-- each line can have different width
-- allowed delimiters: single/multiple spaces, tabs and combination of spaces and tabs


-- Permeability values (md)
-- Simple system 10x1
-- "n*v" represents repeated values : "n" repeated value of "v"

462.7965311791255 2*549.4547649098886

-- comments between lines are allowed

575.4880710664381 469.0640845765829 624.687484573518
741.1518277078972       457.5232943356377
608.3175306515967           
429.1722270419058

/

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
Input file for dead oil PVT table

- Data type

Type    | Size
:---:   | :---:
String | -

```jldoctest
julia> options["PVDO"] = "PVDO.DAT";

```
Example file format:

```

PVDO

-- lines starting with "--" are comments
-- empty lines before keyword are allowed
-- empty lines after / are allowed
-- empty lines between rows are allowed
-- allowed delimiters: single/multiple spaces, tabs and combination of spaces and tabs
-- pref(psi) bo  vis_o(cp)

2465.646 0.9884 1.1663
-- comments between rows are allowed

3335.874 0.9766    1.1695
4138.369254 0.9658 1.1724
4819.467702 0.9565 1.1749    
5801.374962 0.9432 1.1784
6788.358552 0.9299 1.182         
7082.20554 0.9259 1.1831

/

```

---
### PVTW
Input file for the water PVT function

- Data type

Type    | Size
:---:   | :---:
String | -

```jldoctest
julia> options["PVTW"] = "PVTW.DAT";

```

Example file format

```

PVTW

-- lines starting with "--" are comments
-- empty lines before keyword are allowed
-- empty lines after / are allowed
-- empty lines between rows are allowed
-- allowed delimiters: single/multiple spaces, tabs and combination of spaces and tabs
-- see keyword PVTW in Eclipse reference manual

-- pref(psi) bw  cw (1/psi) vis_w(cp) c_visw(1/psi)
3962.194496   1.029     3.17E-06 0.31 0
/


```

---
### SWOF
Input file for the oil-water two-phase relative permeabilty table

- Data type

Type    | Size
:---:   | :---:
String | -

```jldoctest
julia> options["SWOF"] = "SWOF.DAT";

```

Example file format:

```

SWOF

-- lines starting with "--" are comments
-- empty lines before keyword are allowed
-- empty lines after / are allowed
-- empty lines between rows are allowed
-- allowed delimiters: single/multiple spaces, tabs and combination of spaces and tabs

-- sw krw kro pcw(psi)
0.0000 0 0.9000 0
0.1000 0 0.9000 0
0.2000 0.0102 0.6612 0
0.3000 0.0408 0.4592 0
0.4000 0.0918 0.2939 0
0.5000 0.1633 0.1653 0
0.6000 0.2551 0.0735 0
0.7000 0.3673 0.0184 0
0.8000 0.5000 0.0000 0
1.0000 0.5000 0.0000 0.0000
/


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

## Linear solver
---
### linear_solver_backend

- Data type

Type    | Size
:---:   | :---:
String | -

- Supported values
    - `"JULIA"` - linear solvers from [IterativeSolvers.jl](https://github.com/JuliaMath/IterativeSolvers.jl)
    - `"DUNEISTL"` (default) - linear solvers from the [ISTL](https://dune-project.org/modules/dune-istl/) module from [DUNE](https://dune-project.org/)

- Example

```jldoctest
julia> options["linear_solver_backend"] = "JULIA";

julia> options["linear_solver_backend"] = "DUNEISTL";

```

### linear_solver_type

- Data type

Type    | Size
:---:   | :---:
String | -

- Supported values
    - when `options["linear_solver_backend"]="DUNEISTL"`
        - `"GMRES_ILU"` - GMRes solver with incomplete LU preconditioner
        - `"BICGSTAB_ILU"` - BiCGstab solver with incomplete LU preconditioner
    - when `options["linear_solver_backend"]="JULIA"`
        - `"GMRES_ILU"` - GMRes solver with incomplete LU preconditioner
        - `"GMRES_CPR"` - GMRes solver with CPR preconditioner
        - `"BICGSTAB_ILU"` - BiCGstab solver with incomplete LU preconditioner
        - `"BICGSTAB_CPR"` - BiCGstab solver with CPR preconditioner
        - `"JULIA_BACKSLASH"` - Julia backslash solver

- Recommendation for choosing solvers
    - Overall best solver is `"BICGSTAB_ILU"` with `"DUNEISTL"` backend
    - If `"BICGSTAB_ILU"` with `"DUNEISTL"` backend failed
        - Choose `"GMRES_ILU"` or `"BICGSTAB_ILU"` with `"JULIA"` backend for small problems
        - Choose `"GMRES_CPR"` or `"BICGSTAB_CPR"` with `"JULIA"` backend for large problems
        - Only use `"JULIA_BACKSLASH"` for very small problems or problems where robustness is the top priority and other solvers can't converge

- Example

```jldoctest
julia> options["linear_solver_type"] = "BICGSTAB_ILU";

```