using ResSimAD

@doc raw"""
    example1()

Example model one.

# Model description
- Dimension: $30\times15\times1$
- Grid size: $50.0\times50.0\times20.0$ (ft)
- Permeability: 200.0 md
- Porosity: 0.2
- Fluid: Oil-Water
- Initial pressure: 6000.0 psi
- Initial water saturation: 0.1
- Wells:
    - P1: producer at (1, 1, 1), constant BHP 5500.0 psi
    - I1: injector at (30, 15, 1), constant BHP 6500.0 psi
- Schedule:
    - dt0: 0.1 Day
    - dt_max: 50.0 Day
    - t_end: 1825 Day
"""
function example1()
    sim_dir = @__DIR__;
    ## Specify input
    # Grid and Rock
    options = Dict();
    options["nx"] = 30; options["ny"] = 15; options["nz"] = 1;
    options["dx"] = 50.; options["dy"] = 50.; options["dz"] = 20.;
    options["d"] = 8000.;
    options["perm"] = 200.; options["poro"] = 0.2;
    # Fluid
    options["fluid"] = "OW"
    options["po"] = 6000.;
    options["sw"] = 0.1;
    options["PVDO"] = joinpath(sim_dir, "PVDO.DAT");
    options["PVTW"] = joinpath(sim_dir, "PVTW.DAT");
    options["SWOF"] = joinpath(sim_dir, "SWOF.DAT");
    # Wells
    options["producers"] = [];
    p1 = Dict();
    p1["name"] = "P1"; p1["perforation"] = [(1,1,1)]; p1["radius"] = 0.5;
    p1["mode"] = "bhp"; p1["target"] = 5500.;
    push!(options["producers"], p1);

    options["injectors"] = [];
    i1 = Dict();
    i1["name"] = "I1"; i1["perforation"] = [(30,15,1)]; i1["radius"] = 0.5;
    i1["mode"] = "bhp"; i1["target"] = 6500.;
    push!(options["injectors"], i1);
    # Schedule
    options["dt0"] = 0.1
    options["dt_max"] = 50.; options["t_end"] = 10 * 182.5;
    options["min_err"] = 1.0e-3

    sim = Sim(options)
    return sim, options
end
