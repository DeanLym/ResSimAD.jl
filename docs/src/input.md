```@meta
DocTestSetup = quote
    options = Dict()
    Nc = 2
end
```

### Reservoir

##### dx, dy, dz
Cell size along $x$, $y$, $z$ dimension
- Data type

Type    | Size
:---:   | :---:
Float64 | -
Vector{Float64} | $N_c$

- Example:

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


#### Cartesian Grid

##### nx, ny, nz
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
