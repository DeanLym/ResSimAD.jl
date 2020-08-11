

@doc raw"""
    example3()

Example model 3.

# Model description
- Dimension: $60\times60\times1$
- Grid size: $50.0\times50.0\times20.0$ (ft)
- Permeability: Gaussian variogram model
- Porosity: 0.2
- Fluid: Oil-Water
- Initial pressure: 6000.0 psi
- Initial water saturation: 0.1
- Wells:
    - P1: producer at (9, 53, 1), constant BHP 5900.0 psi
    - P2: producer at (51, 11, 1), constant oil rate 500.0 STB/Day
    - P3: producer at (51, 51, 1), constant liquid rate 800.0 STB/Day
    - I1: injector at (8, 24, 1), constant BHP 6100.0 psi
    - I2: injector at (36, 11, 1), constant water rate 800.0 STB/Day
- Schedule:
    - dt0: 0.1 Day
    - dt_max: 50.0 Day
    - t_end: 1825.0 Day
"""
function example3()
    ## Specify input
    # Grid and Rock
    options = Dict();
    options["nx"] = 60; options["ny"] = 60; options["nz"] = 1;
    options["dx"] = 50.; options["dy"] = 50.; options["dz"] = 20.;
    options["d"] = 8000.;
    options["poro"] = 0.2;
    options["perm"] = get_example_data("PERM_GAUSSIAN.DAT");
    # Fluid
    options["fluid"] = "OW"
    options["po"] = 6000.;
    options["sw"] = 0.1;
    options["PVDO"] = get_example_data("PVDO.DAT");
    options["PVTW"] = get_example_data("PVTW.DAT");
    options["SWOF"] = get_example_data("SWOF.DAT");
    # Wells
    options["producers"] = [];
    num_prod = 3
    well_locs = [(9, 53, 1), (51, 11, 1), (51, 51, 1)]
    well_modes = ["bhp", "orat", "lrat"]
    well_targets = [5900., 500., 800.]
    for ind = 1:num_prod
        well = Dict();
        well["name"] = "P"*string(ind);
        well["perforation"] = [well_locs[ind]];
        well["radius"] = 0.5;
        well["mode"] = well_modes[ind];
        well["target"] = well_targets[ind];
        push!(options["producers"], well);
    end

    options["injectors"] = [];
    num_inj = 2
    well_locs = [(8, 24, 1), (36, 11, 1)]
    well_modes = ["bhp", "wrat"]
    well_targets = [6100., -800.]
    for ind = 1:num_inj
        well = Dict();
        well["name"] = "I"*string(ind);
        well["perforation"] = [well_locs[ind]];
        well["radius"] = 0.5;
        well["mode"] = well_modes[ind];
        well["target"] = well_targets[ind];
        push!(options["injectors"], well);
    end

    # Schedule
    options["dt0"] = 0.1
    options["dt_max"] = 50.; options["t_end"] = 10 * 182.5;
    options["min_err"] = 1.0e-3

    sim = Sim(options)
    return sim, options
end
