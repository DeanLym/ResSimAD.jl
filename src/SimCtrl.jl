module SimCtrl

using LinearAlgebra:norm

#! format: off
using ..Global: M, α
using ..AutoDiff: param, data, ones_tensor
using ..Grid: CartGrid, AbstractGrid,
        get_grid_index, set_cell_size, set_perm, set_poro, construct_conn

using ..State: OWState, AbstractState, set_init_state, change_back_state,
        compute_params, update_old_state

using ..Well: StandardWell, WellType, Limit, PRODUCER, INJECTOR,
        get_ctrl_mode, compute_wi, compute_qo, compute_qw

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
    function Sim(nx::Int, ny::Int, nz::Int)::Sim
        grid = CartGrid(nx, ny, nz)
        state = OWState(grid.numcell, grid.connlist.numconn)
        scheduler = Scheduler()
        nsolver = NRSolver()
        injectors = Vector{StandardWell{INJECTOR}}()
        producers = Vector{StandardWell{PRODUCER}}()
        return new(grid, state, scheduler, nsolver, injectors, producers)
    end
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

function setup(sim::Sim, input::Dict)::Nothing
    grid, state, scheduler = sim.grid, sim.state, sim.scheduler
    producers, injectors = sim.producers, sim.injectors
    # Set up grid
    dx, dy, dz = input["dx"], input["dy"], input["dz"]
    set_cell_size(grid, dx, dy, dz)
    set_perm(grid, input["k"])
    set_poro(grid, input["ϕ"])
    construct_conn(grid)

    # Set up initial state
    set_init_state(state, input["p0"], input["sw0"])

    # Set up wells
    numvar = state.numvar
    for p in input["producers"]
        push!(producers, init_well(PRODUCER, p, numvar, grid))
    end

    for p in input["injectors"]
        push!(injectors, init_well(INJECTOR, p, numvar, grid))
    end

    # Set up options
    scheduler.t_current = 0.0
    scheduler.t0 = 0.0
    scheduler.dt_max = get(input, "dt_max", 100.)
    scheduler.t_end = get(input, "t_end", 100.)

    return nothing
end


function compute_well_rate(sim::Sim)
    producers, injectors = sim.producers, sim.injectors
    state = sim.state
    for prod in producers
        compute_qo(prod, state)
        compute_qw(prod, state)
    end
    for inj in injectors
        compute_qw(inj, state)
    end
end


function runsim(sim::Sim)
    grid, state, sch, nsolver = sim.grid, sim.state, sim.scheduler, sim.nsolver
    producers, injectors = sim.producers, sim.injectors
    p_all = Vector{Vector{Float64}}()
    sw_all = Vector{Vector{Float64}}()

    t_res, t_jac, t_sol = 0.0, 0.0, 0.0
    tmp = nothing
    num_newton = 0
    while sch.t_current < sch.t_end
        println(sch.t_current)
        # println("!=================== Timestep $it ====================!")
        # Newton Iteration
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
        if converge
            update_old_state(state)
        else
            change_back_state(state)
            compute_params(state)
        end
        # Save solution
        push!(p_all, data(state.p))
        push!(sw_all, data(state.sw))
    end
    println("\nt_res： $t_res, t_jac: $t_jac, t_sol: $t_sol")
    println("NumNewton: $num_newton\n")
    return p_all, sw_all
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
