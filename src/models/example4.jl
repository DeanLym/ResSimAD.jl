

@doc raw"""
    example4()

Example model 4.

# Model description
- Dimension: $60\times60\times40$
- Grid size: $50.0\times50.0\times20.0$ (ft)
- Permeability:
    - kx: 2000.0 md
    - ky: 2000.0 md
    - kz: 500.0 md
- Porosity: 0.2
- Fluid: Oil-Water
- Initial pressure: 6000.0 psi at top layer
- Initial water saturation: 0.1
- Wells:
    - P1: producer at (9, 53, 1), constant BHP 5900.0 psi
    - P2: producer at (51, 11, 1), constant oil rate 500.0 STB/Day
    - I1: injector at (8, 24, 1), constant BHP 6100.0 psi
    - I2: injector at (36, 11, 1), constant water rate 800.0 STB/Day
- Schedule:
    - dt0: 0.01 Day
    - dt_max: 30.0 Day
    - t_end: 1825.0 Day
"""
function example4()
    ## Specify input
    # Grid and Rock
    options = Dict();
    options["nx"] = 60; options["ny"] = 60; options["nz"] = 40;
    options["dx"] = 50.; options["dy"] = 50.; options["dz"] = 20.;
    options["tops"] = 8000.;
    options["d"] = 10000.0;
    options["poro"] = 0.2;
    options["perm"] = 2000.0;
    options["multpermz"] = 0.25;
    # Fluid
    options["fluid"] = "OW"
    options["equil"] = (8000.0, 6000.0);
    # options["po"] = 6000.0;
    options["sw"] = 0.1;
    options["PVDO"] = get_example_data("PVDO.DAT");
    options["PVTW"] = get_example_data("PVTW.DAT");
    options["SWOF"] = get_example_data("SWOF.DAT");
    # Wells
    options["producers"] = [];
    num_prod = 2
    well_locs = [(13, 7, 10), (13, 21, 10)]
    for ind = 1:num_prod
        well = Dict();
        well["name"] = "P"*string(ind);
        well["perforation"] = [well_locs[ind]];
        well["radius"] = 0.5;
        well["mode"] = "bhp";
        well["target"] = 5400.0;
        push!(options["producers"], well);
    end

    options["injectors"] = [];
    num_inj = 2
    well_locs = [(25, 21, 30), (25, 53, 30)]
    for ind = 1:num_inj
        well = Dict();
        well["name"] = "I"*string(ind);
        well["perforation"] = [well_locs[ind]];
        well["radius"] = 0.5;
        well["mode"] = "bhp";
        well["target"] = 6600.0;
        push!(options["injectors"], well);
    end

    # Nonlinear solver options
    options["max_newton_iter"] = 15
    options["min_err"] = 1.0e-3

    # Linear solver options
    options["linear_solver"] = "GMRES_CPR"

    # Schedule
    options["dt0"] = 0.01
    options["dt_max"] = 15.; options["t_end"] = 800.0;

    sim = Sim(options)
    return sim, options
end
