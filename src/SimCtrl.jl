module SimCtrl

using LinearAlgebra:norm
using IterativeSolvers, IncompleteLU

#! format: off
using ..Global: M, α, β, g_, gc
using ..AutoDiff: param, data, ones_tensor
using ..Grid: CartGrid, AbstractGrid, ConnList, set_cell_depth,
        get_grid_index, set_cell_size, set_perm, set_poro, construct_connlist

using ..Fluid: AbstractPVTO, AbstractPVTW, AbstractSWOF, SWOFTable, PVDO, PVTW

using ..State: OWState, AbstractState, set_init_state, change_back_state,
        compute_params, update_old_state, save_state_result

using ..Well: StandardWell, WellType, Limit, PRODUCER, INJECTOR,
        get_ctrl_mode, compute_wi, compute_qo, compute_qw, save_result

using ..Schedule: Scheduler, update_dt, reset_dt, insert_time_step, set_time_step

using ..Solver: NonlinearSolver, NRSolver, compute_residual_error,
        assemble_residual, assemble_jacobian, update_solution
#! format: on

struct Sim
    grid::CartGrid
    state::OWState
    scheduler::Scheduler
    nsolver::NonlinearSolver
    injectors::Dict{String, StandardWell{INJECTOR}}
    producers::Dict{String, StandardWell{PRODUCER}}
end

get_well_type =
    Dict{String,WellType}("producer" => PRODUCER, "injector" => INJECTOR)

function Sim(input::Dict)::Sim
    # Construct Grid
    nx, ny, nz = input["nx"], input["ny"], input["nz"]
    grid = CartGrid(nx, ny, nz)
    dx, dy, dz = input["dx"], input["dy"], input["dz"]
    set_cell_size(grid, dx, dy, dz)
    set_cell_depth(grid, input["d"])
    set_perm(grid, input["k"])
    set_poro(grid, input["ϕ"])
    construct_connlist(grid)

    # Construct Fluid
    if "SWOF" in keys(input)
        swof = SWOFTable(input["SWOF"])
    elseif "SWOFCorey" in keys(input)
        swof = SWOFCorey(input["SWOFCorey"])
    end
    pvdo = PVDO(input["PVDO"])
    pvtw = PVTW(input["PVTW"])

    # Construct State
    ρo_std = get(input, "ρo_std", 49.1)
    ρw_std = get(input, "ρw_std", 64.79)
    #! format: off
    state = OWState(grid.numcell, grid.connlist.numconn,
                      ρo_std, ρw_std, pvdo, pvtw, swof)
    #! format: on
    set_init_state(state, input["p0"], input["sw0"])

    # Construct Wells
    producers = Dict{String, StandardWell{PRODUCER}}()
    injectors = Dict{String, StandardWell{INJECTOR}}()
    well_names = Set{String}()
    numvar = state.numvar
    for p in input["producers"]
        name = p["name"]
        @assert !(name in well_names)
        push!(well_names, name)
        producers[name] = init_well(PRODUCER, p, numvar, grid)
    end

    for p in input["injectors"]
        name = p["name"]
        @assert !(name in well_names)
        push!(well_names, name)
        injectors[p["name"]] = init_well(INJECTOR, p, numvar, grid)
    end

    # Construct Scheduler
    sch = Scheduler()
    sch.t_current = 0.0
    sch.t0 = 0.0
    sch.dt_max = get(input, "dt_max", 100.)
    sch.dt0 = get(input, "dt0", 0.01)
    reset_dt(sch)
    t_end = get(input, "t_end", 100.)
    set_time_step(sch, get(input, "time_step", [t_end]))

    # Construct Nonlinear Solver
    nsolver = NRSolver()
    nsolver.max_newton_iter = get(input, "max_newton_iter", 10)
    nsolver.min_err = get(input, "min_err", 1.0e-6)

    return Sim(grid, state, sch, nsolver, injectors, producers)
end


function init_well(T::WellType, well_option::Dict, numvar::Int, grid::CartGrid)
    name = well_option["name"]
    perf = Vector{Int}()
    for indices in well_option["perforation"]
        push!(perf, get_grid_index(grid, indices...))
    end
    radius = get(well_option, "radius", 0.5)
    #
    well = StandardWell{T}(name, perf, radius, numvar)
    # Set control mode and target
    well.mode = get_ctrl_mode[well_option["mode"]]
    well.target = well_option["target"]
    # Set limits
    limits = get(well_option, "limits", Dict{Limit, Float64}())
    for (k,v) in limits
        well.limits[k] = v
    end
    # Comput well index
    compute_wi(well, grid)
    return well
end

function add_well(sim::Sim, welltype::String, well_option::Dict)::Nothing
    name = well_option["name"]
    @assert !(name in keys(sim.producers)) && !(name in keys(sim.injectors))
    T = get_well_type[lowercase(welltype)]
    nv, grid = sim.state.numvar, sim.grid
    if T == PRODUCER
        sim.producers[name] = init_well(T, well_option, nv, grid)
    elseif T == INJECTOR
        sim.injectors[name] = init_well(T, well_option, nv, grid)
    end
    reset_dt(sim.scheduler)
    return nothing
end

function change_well_mode(sim::Sim, name::String, mode::String, target::Float64)::Nothing
    @assert name in keys(sim.producers) || name in keys(sim.injectors)
    if name in keys(sim.producers)
        sim.producers[name].mode = get_ctrl_mode[mode]
        sim.producers[name].target = target
    else
        sim.injectors[name].mode = get_ctrl_mode[mode]
        sim.injectors[name].target = target
    end
    reset_dt(sim.scheduler)
    return nothing
end

function change_well_target(sim::Sim, name::String, target::Float64; cut_dt=false)::Nothing
    @assert name in keys(sim.producers) || name in keys(sim.injectors)
    if name in keys(sim.producers)
        sim.producers[name].target = target
    else
        sim.injectors[name].target = target
    end
    if cut_dt
        reset_dt(sim.scheduler)
    end
    return nothing
end

function get_well_rates(sim::Sim, name::String, data::String)
    @assert name in keys(sim.producers) || name in keys(sim.injectors)
    if name in keys(sim.producers)
        return sim.producers[name].results[!, data]
    else
        return sim.injectors[name].results[!, data]
    end
end

function shut_well(sim::Sim, name::String)::Nothing
    @assert name in keys(sim.producers) || name in keys(sim.injectors)
    if name in keys(sim.producers)
        sim.producers[name].mode = get_ctrl_mode["shut"]
    else
        sim.injectors[name].mode = get_ctrl_mode["shut"]
    end
    return nothing
end


function compute_well_rate(sim::Sim)
    state = sim.state
    for prod in values(sim.producers)
        compute_qo(prod, state)
        compute_qw(prod, state)
    end
    for inj in values(sim.injectors)
        compute_qw(inj, state)
    end
end

function save_well_results(sim::Sim, t::Float64)
    for prod in values(sim.producers)
        save_result(prod, t)
    end
    for inj in values(sim.injectors)
        save_result(inj, t)
    end
end

function newton_step(sim::Sim)::Nothing
    grid, state, sch, nsolver = sim.grid, sim.state, sim.scheduler, sim.nsolver
    producers, injectors = sim.producers, sim.injectors
    dt = sch.dt
    compute_well_rate(sim)
    compute_ro(state, grid, producers, dt)
    compute_rw(state, grid, producers, injectors, dt)
    err = compute_residual_error(state, grid, dt)
    if err < nsolver.min_err
        nsolver.converged = true
        return nothing
    end
    nsolver.residual = assemble_residual(state)
    # Compute jacobian
    nsolver.jac = assemble_jacobian(state)
    # Solve equation
    # nsolver.δx = nsolver.jac \ nsolver.residual
    nsolver.δx = gmres(nsolver.jac, nsolver.residual, Pl=ilu(nsolver.jac, τ=0.1))
    update_solution(state, nsolver.δx)
    compute_params(state)
    nsolver.converged = false
    return nothing
end


function step(sim::Sim)::Nothing
    grid, state, sch, nsolver = sim.grid, sim.state, sim.scheduler, sim.nsolver
    println(sch.t_next)
    nsolver.converged, err, newton_iter = false, Inf, 0

    while !nsolver.converged
        if newton_iter > nsolver.max_newton_iter
            break
        end
        newton_step(sim)
        newton_iter += 1
    end
    update_dt(sch, state, nsolver.converged)
    push!(nsolver.num_iter, newton_iter)
    if nsolver.converged
        save_state_result(state, sch.t_current)
        update_old_state(state)
        save_well_results(sim, sch.t_current)
    else
        change_back_state(state)
        compute_params(state)
    end
    println("NumNewton: $newton_iter\n")
    return nothing
end

function step_to(sim::Sim, t::Float64)::Nothing
    sch = sim.scheduler
    insert_time_step(sch, t)
    while sch.t_current < t
        step(sim)
    end
end

function runsim(sim::Sim)::Nothing
    sch = sim.scheduler
    while sch.t_current < sch.time_step[end]
        step(sim)
    end
    return nothing
end

#! format: off
function compute_ao(state::OWState, grid::AbstractGrid, dt::Float64)::Nothing
    state.ao .= 1.0 / M .* grid.v .* grid.ϕ .*
        (state.so ./ state.bo - state.son ./ state.bon) / dt
    return nothing
end
#! format: on

#! format: off
function compute_aw(state::OWState, grid::AbstractGrid, dt::Float64)::Nothing
    state.aw .= 1.0 / M .* grid.v .* grid.ϕ .*
        (state.sw ./ state.bw - state.swn ./ state.bwn) / dt
    return nothing
end
#! format: on

function compute_γo(state::OWState, connlist::ConnList)::Nothing
    l_list, r_list = connlist.l, connlist.r
    ρo, γo = state.ρo, state.γo
    for i = 1:connlist.numconn
        l, r = l_list[i], r_list[i]
        γo[i] = β * g_ / gc * (ρo[l] + ρo[r]) / 2.
    end
    return nothing
end

function compute_γw(state::OWState, connlist::ConnList)::Nothing
    l_list, r_list = connlist.l, connlist.r
    ρw, γw = state.ρw, state.γw
    for i = 1:connlist.numconn
        l, r = l_list[i], r_list[i]
        γw[i] = β * g_ / gc * (ρw[l] + ρw[r]) / 2.
    end
    return nothing
end

function compute_ΔΨo(state::OWState, connlist::ConnList)::Nothing
    l_list, r_list, Δd = connlist.l, connlist.r, connlist.Δd
    ΔΨo, p, γo = state.ΔΨo, state.p, state.γo
    for i = 1:connlist.numconn
        l, r = l_list[i], r_list[i]
        ΔΨo[i] = p[l] - p[r] - γo[i]*Δd[i]
    end
    return nothing
end

function compute_ΔΨw(state::OWState, connlist::ConnList)::Nothing
    l_list, r_list, Δd = connlist.l, connlist.r, connlist.Δd
    ΔΨw, p, γw = state.ΔΨw, state.p, state.γw
    for i = 1:connlist.numconn
        l, r = l_list[i], r_list[i]
        ΔΨw[i] = p[l] - p[r] - γw[i]*Δd[i]
    end
    return nothing
end

function compute_fo(state::OWState, connlist::ConnList)::Nothing
    l_list, r_list, trans = connlist.l, connlist.r, connlist.trans
    ΔΨo, λo, fo = state.ΔΨo, state.λo, state.fo
    for i = 1:connlist.numconn
        l, r = l_list[i], r_list[i]
        if ΔΨo[i] > 0.0
            fo[i] = λo[l] * trans[i] * ΔΨo[i]
        else
            fo[i] = λo[r] * trans[i] * ΔΨo[i]
        end
    end
    return nothing
end

function compute_fw(state::OWState, connlist::ConnList)::Nothing
    l_list, r_list, trans = connlist.l, connlist.r, connlist.trans
    ΔΨw, λw, fw = state.ΔΨw, state.λw, state.fw
    for i = 1:connlist.numconn
        l, r = l_list[i], r_list[i]
        if ΔΨw[i] > 0.0
            fw[i] = λw[l] * trans[i] * ΔΨw[i]
        else
            fw[i] = λw[r] * trans[i] * ΔΨw[i]
        end
    end
    return nothing
end


function compute_ro(
    state::OWState,
    grid::AbstractGrid,
    producers::Dict{String, StandardWell{PRODUCER}},
    dt::Float64,
)::Nothing
    connlist = grid.connlist
    compute_ao(state, grid, dt)
    compute_γo(state, connlist)
    compute_ΔΨo(state, connlist)
    compute_fo(state, connlist)
    # Accumulation term
    state.ro .= - state.ao
    # Sink / source term
    for w in values(producers)
        state.ro[w.ind] -= w.qo
    end
    # Flux term
    l, r = connlist.l, connlist.r
    for i = 1:connlist.numconn
        state.ro[l[i]] -= state.fo[i]
        state.ro[r[i]] += state.fo[i]
    end
    return nothing
end

function compute_rw(
    state::OWState,
    grid::AbstractGrid,
    producers::Dict{String, StandardWell{PRODUCER}},
    injectors::Dict{String, StandardWell{INJECTOR}},
    dt::Float64,
)::Nothing
    connlist = grid.connlist
    compute_aw(state, grid, dt)
    compute_γw(state, connlist)
    compute_ΔΨw(state, connlist)
    compute_fw(state, connlist)
    # Accumulation term
    state.rw .= - state.aw
    # Sink / source term
    for w in values(producers)
        state.rw[w.ind] -= w.qw
    end
    for w in values(injectors)
        state.rw[w.ind] -= w.qw
    end
    # Flux term
    l, r = connlist.l, connlist.r
    for i = 1:grid.connlist.numconn
        state.rw[l[i]] -= state.fw[i]
        state.rw[r[i]] += state.fw[i]
    end
    return nothing
end


end
