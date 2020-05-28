module SimCtrl

using LinearAlgebra:norm

#! format: off
using ..AutoDiff: param, data
using ..Grid: CartGrid, AbstractGrid,
        get_grid_index, set_cell_size, set_perm, set_poro, construct_conn

using ..State: OWState, AbstractState, set_init_state,
        compute_ro, compute_rw, compute_params, update_old_state

using ..Schedule: Scheduler, update_dt

using ..Solver: NonlinearSolver, NRSolver, compute_residual_error,
        assemble_residual, assemble_jacobian, update_solution
#! format: on

struct Sim
    grid::CartGrid
    state::OWState
    scheduler::Scheduler
    nsolver::NonlinearSolver
    inj_bhp::Vector{Float64}
    prod_bhp::Vector{Float64}
    function Sim(nx::Int, ny::Int, nz::Int)::Sim
        grid = CartGrid(nx, ny, nz)
        state = OWState(grid.numcell, grid.connlist.numconn)
        scheduler = Scheduler()
        nsolver = NRSolver()
        inj_bhp = Inf * ones(grid.numcell)
        prod_bhp = Inf * ones(grid.numcell)
        return new(grid, state, scheduler, nsolver, inj_bhp, prod_bhp)
    end
end

function setup(sim::Sim, input::Dict{Any, Any})::Nothing
    # Set up grid
    dx, dy, dz = input["dx"], input["dy"], input["dz"]
    set_cell_size(sim.grid, dx, dy, dz)
    set_perm(sim.grid, input["k"])
    set_poro(sim.grid, input["ϕ"])
    construct_conn(sim.grid)

    # Set up initial state
    set_init_state(sim.state, input["p0"], input["sw0"])

    # Set up wells
    for (x, y, z, bhp) in input["producer"]
        ind = get_grid_index(sim.grid, x, y, z)
        sim.prod_bhp[ind] = bhp
    end
    for (x, y, z, bhp) in input["injector"]
        ind = get_grid_index(sim.grid, x, y, z)
        sim.inj_bhp[ind] = bhp
    end

    return nothing
end


function runsim(sim::Sim)
    grid, state, sch, nsolver = sim.grid, sim.state, sim.scheduler, sim.nsolver
    prod_bhp, inj_bhp = sim.prod_bhp, sim.inj_bhp
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
            compute_ro(state, grid, prod_bhp, dt)
            compute_rw(state, grid, prod_bhp, inj_bhp, dt)
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



end
