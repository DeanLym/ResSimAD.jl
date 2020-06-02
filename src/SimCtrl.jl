module SimCtrl

using LinearAlgebra:norm

#! format: off
using ..Global: M, α
using ..AutoDiff: param, data, ones_tensor
using ..Grid: CartGrid, AbstractGrid,
        get_grid_index, set_cell_size, set_perm, set_poro, construct_connlist

using ..Fluid: AbstractPVTO, AbstractPVTW, AbstractSWOF, SWOFTable, PVDO, PVTW

using ..State: OWState, AbstractState, set_init_state, change_back_state,
        compute_params, update_old_state, save_state_result

using ..Well: StandardWell, WellType, Limit, PRODUCER, INJECTOR,
        get_ctrl_mode, compute_wi, compute_qo, compute_qw, save_result

using ..Schedule: Scheduler, update_dt

using ..Solver: NonlinearSolver, NRSolver, compute_residual_error,
        assemble_residual, assemble_jacobian, update_solution
#! format: on

struct Sim
    grid::CartGrid
    state::OWState
    scheduler::Scheduler
    nsolver::NonlinearSolver
    injectors::Vector{StandardWell{INJECTOR}}
    producers::Vector{StandardWell{PRODUCER}}
end

function Sim(input::Dict)::Sim
    # Construct Grid
    nx, ny, nz = input["nx"], input["ny"], input["nz"]
    grid = CartGrid(nx, ny, nz)
    dx, dy, dz = input["dx"], input["dy"], input["dz"]
    set_cell_size(grid, dx, dy, dz)
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
    state = OWState(grid.numcell, grid.connlist.numconn, pvdo, pvtw, swof)
    set_init_state(state, input["p0"], input["sw0"])

    # Construct Wells
    injectors = Vector{StandardWell{INJECTOR}}()
    producers = Vector{StandardWell{PRODUCER}}()
    numvar = state.numvar
    for p in input["producers"]
        push!(producers, init_well(PRODUCER, p, numvar, grid))
    end

    for p in input["injectors"]
        push!(injectors, init_well(INJECTOR, p, numvar, grid))
    end

    # Construct Scheduler
    scheduler = Scheduler()
    scheduler.t_current = 0.0
    scheduler.t0 = 0.0
    scheduler.dt_max = get(input, "dt_max", 100.)
    scheduler.t_end = get(input, "t_end", 100.)
    scheduler.report_time = get(input, "report_time", [scheduler.t_end])

    # Construct Nonlinear Solver
    nsolver = NRSolver()

    return Sim(grid, state, scheduler, nsolver, injectors, producers)
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

function compute_well_rate(sim::Sim)
    state = sim.state
    for prod in sim.producers
        compute_qo(prod, state)
        compute_qw(prod, state)
    end
    for inj in sim.injectors
        compute_qw(inj, state)
    end
end

function save_well_results(sim::Sim, t::Float64)
    for prod in sim.producers
        save_result(prod, t)
    end
    for inj in sim.injectors
        save_result(inj, t)
    end
end

function step(sim::Sim)::Nothing
    grid, state, sch, nsolver = sim.grid, sim.state, sim.scheduler, sim.nsolver
    producers, injectors = sim.producers, sim.injectors
    println(sch.t_next)
    t_res, t_jac, t_sol = 0.0, 0.0, 0.0
    num_newton = 0
    dt = sch.dt
    converge, err, newton_iter = true, Inf, 0
    while err > nsolver.min_err
        #print("Netwon Iteration $newton_iter: ")
        # newton_iter += 1
        # Calculate residual and jacobian
        t0 = time()
        compute_well_rate(sim)
        compute_ro(state, grid, producers, dt)
        compute_rw(state, grid, producers, injectors, dt)
        err = compute_residual_error(state, grid, dt)
        if err < nsolver.min_err
            num_newton += newton_iter
            break
        end
        residual = assemble_residual(state)
        tt = time() - t0
        #print(" t_res: ", tt)
        t_res += tt

        # Compute jacobian
        t0 = time()
        jac = assemble_jacobian(state)
        tt = time() - t0
        #print(" t_jac: ", tt)
        t_jac += tt

        # Solve equation
        t0 = time()
        δx = jac \ residual
        tt = time() - t0
        #print(" t_sol: $tt \n")
        t_sol += tt
        newton_iter += 1
        # Update variables
        if newton_iter > nsolver.max_iter
            converge = false
            break
        else
            update_solution(state, δx)
            compute_params(state)
        end
    end
    update_dt(sch, state, converge)
    push!(nsolver.num_iter, num_newton)
    if converge
        save_state_result(state, sch.t_next)
        update_old_state(state)
        save_well_results(sim, sch.t_next)
    else
        change_back_state(state)
        compute_params(state)
    end
    println("t_res： $t_res, t_jac: $t_jac, t_sol: $t_sol")
    println("NumNewton: $num_newton\n")
    return nothing
end

function runsim(sim::Sim)::Nothing
    sch = sim.scheduler
    while sch.t_current < sch.t_end
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

function compute_fo(state::OWState, grid::AbstractGrid)::Nothing
    trans = grid.connlist.trans
    for i = 1:grid.connlist.numconn
        l, r = grid.connlist.l[i], grid.connlist.r[i]
        if state.p[l] > state.p[r]
            state.fo[i] = state.λo[l] * trans[i] * (state.p[l] - state.p[r])
        else
            state.fo[i] = state.λo[r] * trans[i] * (state.p[l] - state.p[r])
        end
    end
    return nothing
end

function compute_fw(state::OWState, grid::AbstractGrid)::Nothing
    trans = grid.connlist.trans
    for i = 1:grid.connlist.numconn
        l, r = grid.connlist.l[i], grid.connlist.r[i]
        if state.p[l] > state.p[r]
            state.fw[i] = state.λw[l] * trans[i] * (state.p[l] - state.p[r])
        else
            state.fw[i] = state.λw[r] * trans[i] * (state.p[l] - state.p[r])
        end
    end
    return nothing
end

function compute_ro(
    state::OWState,
    grid::AbstractGrid,
    producers::Vector{StandardWell{PRODUCER}},
    dt::Float64,
)::Nothing
    compute_ao(state, grid, dt)
    compute_fo(state, grid)
    # Accumulation term
    state.ro .= - state.ao
    # Sink / source term
    for w in producers
        state.ro[w.ind] -= w.qo
    end
    # Flux term
    l, r = grid.connlist.l, grid.connlist.r
    for i = 1:grid.connlist.numconn
        state.ro[l[i]] -= state.fo[i]
        state.ro[r[i]] += state.fo[i]
    end
    return nothing
end

function compute_rw(
    state::OWState,
    grid::AbstractGrid,
    producers::Vector{StandardWell{PRODUCER}},
    injectors::Vector{StandardWell{INJECTOR}},
    dt::Float64,
)::Nothing
    compute_aw(state, grid, dt)
    compute_fw(state, grid)
    # Accumulation term
    state.rw .= - state.aw
    # Sink / source term
    for w in producers
        state.rw[w.ind] -= w.qw
    end
    for w in injectors
        state.rw[w.ind] -= w.qw
    end
    # Flux term
    l, r = grid.connlist.l, grid.connlist.r
    for i = 1:grid.connlist.numconn
        state.rw[l[i]] -= state.fw[i]
        state.rw[r[i]] += state.fw[i]
    end
    return nothing
end


end
