module SimCtrl

using LinearAlgebra:norm

#! format: off
using ..AutoDiff: param
using ..Grid: CartGrid, AbstractGrid,
        get_grid_index, set_cell_size, set_perm, set_poro, construct_conn

using ..State: OWState, AbstractState, set_init_state,
        compute_ro, compute_rw, compute_params, update_old_state

using ..Solver:assemble_residual, assemble_jacobian, update_solution
#! format: on

struct Sim
    grid::CartGrid
    state::OWState
    inj_bhp::Vector{Float64}
    prod_bhp::Vector{Float64}
    function Sim(nx::Int, ny::Int, nz::Int)::Sim
        grid = CartGrid(nx, ny, nz)
        state = OWState(grid.numcell, grid.connlist.numconn)
        inj_bhp = Inf * ones(grid.numcell)
        prod_bhp = Inf * ones(grid.numcell)
        return new(grid, state, inj_bhp, prod_bhp)
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


function runsim(sim::Sim, numtime::Int)::Tuple{Matrix{Float64}, Matrix{Float64}}
    dt = 1.0 # Hard code
    p_all = zeros(numtime, sim.grid.numcell)
    sw_all = zeros(numtime, sim.grid.numcell)

    t_res, t_jac, t_sol = 0.0, 0.0, 0.0
    tmp = nothing
    num_newton = 0
    for it=1:numtime
        println("!=================== Timestep $it ====================!")
        # Newton Iteration
        min_err = 1e-6
        err = min_err + 1
        newton_iter = 0
        while err > min_err
            print("Netwon Iteration $newton_iter: ")
            # newton_iter += 1
            # Calculate residual and jacobian
            t0 = time()
            compute_ro(sim.state, sim.grid, sim.prod_bhp, dt)
            compute_rw(sim.state, sim.grid, sim.prod_bhp, sim.inj_bhp, dt)
            residual = assemble_residual(sim.state)
            tt = time() - t0
            print(" t_res: ", tt)
            t_res += tt

            err = norm(residual)
            if err < min_err
                println("")
                num_newton += newton_iter
                break
            end
            # Compute jacobian
            t0 = time()
            jac = assemble_jacobian(sim.state)
            tt = time() - t0
            print(" t_jac: ", tt)
            t_jac += tt

            # Solve equation
            t0 = time()
            δx = jac \ residual
            tt = time() - t0
            print(" t_sol: $tt \n")
            t_sol += tt
            newton_iter += 1
            # Update variables
            update_solution(sim.state, δx)
            compute_params(sim.state)
        end
        # Update state
        update_old_state(sim.state)
        # Save solution
        for ind=1:sim.grid.numcell
            p_all[it, ind] = sim.state.p[ind].val
            sw_all[it, ind] = sim.state.sw[ind].val
        end
    end
    println("\nt_res： $t_res, t_jac: $t_jac, t_sol: $t_sol")
    println("NumNewton: $num_newton\n")
    return p_all, sw_all
end



end
