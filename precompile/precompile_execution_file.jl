using ResSimAD

sim, options = get_model("example1");

change_well_target(sim, "P1", 5800.);
change_well_mode(sim, "I1", "bhp", 6300.);

well = Dict("name" => "I2", "perforation"=>[(8,12,1)],
          "radius"=>0.5, "mode"=>"bhp", "target"=>6500.);
add_well(sim, "injector", well);

shut_well(sim, "I1");

change_dt(sim, 1.0);

step_to(sim, 100.);

runsim(sim);

t = get_well_rates(sim, "P1", "TIME");
po = get_state_map(sim, "po", t[end]);
