

@doc raw"""
    example2()

Example model one.

# Model description
- Dimension: $60\times60\times1$
- Grid size: $200.0\times200.0\times60.0$ (ft)
- Permeability: bimodal channel system
- Porosity: 0.2
- Fluid: Oil-Water
- Initial pressure: 6000.0 psi
- Initial water saturation: 0.1
- Wells:
    - P1: producer at (20, 18, 1), constant BHP 5900.0 psi
    - P2: producer at (13, 35, 1), constant BHP 5900.0 psi
    - I1: injector at (40, 40, 1), constant BHP 6100.0 psi
    - I2: injector at (25, 53, 1), constant BHP 6100.0 psi
- Schedule:
    - dt0: 2.0 Day
    - dt_max: 30.0 Day
    - t_end: 1825.0 Day
"""
function example2()
    ## Specify input
    # Grid and Rock
    options = Dict();
    options["nx"] = 60; options["ny"] = 60; options["nz"] = 1;
    options["dx"] = 200.; options["dy"] = 200.; options["dz"] = 60.;
    options["d"] = 8000.;
    options["poro"] = 0.2;
    options["perm"] = get_example_data("PERM_BIMODAL.DAT");
    # Fluid
    options["fluid"] = "OW"
    options["po"] = 6000.;
    options["sw"] = 0.1;
    options["PVDO"] = get_example_data("PVDO.DAT");
    options["PVTW"] = get_example_data("PVTW.DAT");
    options["SWOF"] = get_example_data("SWOF.DAT");
    # Wells
    options["producers"] = [];
    w = Dict("name"=>"P1", "perforation"=>[(20,18,1)],
              "radius"=>0.5, "mode"=>"bhp", "target"=>5900.);
    push!(options["producers"], w);
    w = Dict("name"=>"P2", "perforation"=>[(13,35,1)],
              "radius"=>0.5, "mode"=>"bhp", "target"=>5900.);
    push!(options["producers"], w);

    options["injectors"] = [];
    w = Dict("name"=>"I1", "perforation"=>[(40,40,1)],
              "radius"=>0.5, "mode"=>"bhp", "target"=>6100.);
    push!(options["injectors"], w);
    w = Dict("name" => "I2", "perforation"=>[(25,53,1)],
              "radius"=>0.5, "mode"=>"bhp", "target"=>6100.);
    push!(options["injectors"], w);

    # Schedule
    options["dt0"] = 2.0
    options["dt_max"] = 30.; options["t_end"] = 10 * 182.5;
    options["min_err"] = 1.0e-3

    sim = Sim(options)
    return sim, options
end
