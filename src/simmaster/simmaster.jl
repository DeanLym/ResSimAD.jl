module SimMaster

using LinearAlgebra:norm
using SparseArrays:sparse, SparseMatrixCSC

#! format: off
using ..Global: M, α, β, g_, gc
using ..AutoDiff: param, data, ones_tensor, zeros_tensor, grad
using ..InputParse: parse_input_grid, parse_input_rock, parse_input_fluid

using ..Grid: CartGrid, AbstractGrid, construct_connlist, set_cell_depth,
                get_grid_index, set_cell_size
using ..Rock: AbstractRock, StandardRock, set_perm, set_poro
using ..Fluid: AbstractFluid, OWFluid, PVT, PVTC, SWOFTable, SWOFCorey,
                set_fluid_tn, update_phases, compute_a, update_primary_variable,
                update_fluid_tn, reset_primary_variable
using ..Reservoir: AbstractReservoir, StandardReservoir, save_fluid_results

using ..Facility: AbstractFacility, StandardWell, WellType, Limit, PRODUCER, INJECTOR,
        get_ctrl_mode, compute_wi, compute_qo, compute_qw, save_result

using ..Schedule: Scheduler, update_dt, reset_dt, insert_time_step, set_time_step

using ..LinearSolver: AbstractLinearSolver, Julia_BackSlash_Solver, GMRES_ILU0_Solver, solve
#! format: on

abstract type NonlinearSolver end

struct Sim
    reservoir::AbstractReservoir
    facility::Dict{String, AbstractFacility}
    scheduler::Scheduler
    nsolver::NonlinearSolver
    lsolver::AbstractLinearSolver
end

# Newton Raphson Solver
mutable struct NRSolver <: NonlinearSolver
    max_newton_iter::Int
    newton_iter::Int
    min_err::Float64
    num_iter::Vector{Int}
    converged::Bool
    δx::Vector{Float64}
    residual::Vector{Float64}
    jac::SparseMatrixCSC{Float64,Int}
    NRSolver() = new(10, 0, 1.0e-6, Vector{Int}(), false)
end

function newton_step(sim::Sim)::Nothing
    reservoir = sim.reservoir
    grid, fluid, rock = reservoir.grid, reservoir.fluid, reservoir.rock
    facility, nsolver, lsolver, sch = sim.facility, sim.nsolver, sim.lsolver, sim.scheduler
    # Assemble residual
    nsolver.residual = assemble_residual(fluid)
    # Assemble jacobian
    nsolver.jac = assemble_jacobian(fluid)
    # Solve equation
    nsolver.δx = solve(lsolver, nsolver.jac, nsolver.residual)
    # Update solution
    update_solution(fluid, nsolver.δx)
    # Update dynamic states and then compute new residual
    update_phases(fluid, grid.connlist)
    compute_residual(fluid, grid, rock, facility, sch.dt)
    nsolver.newton_iter += 1
    return nothing
end

function step(sim::Sim)::Nothing
    reservoir = sim.reservoir
    grid, fluid, rock = reservoir.grid, reservoir.fluid, reservoir.rock
    facility, nsolver, lsolver, sch = sim.facility, sim.nsolver, sim.lsolver, sim.scheduler
    println(sch.t_next)
    while true
        # Check convergence
        err = compute_residual_error(fluid, grid, rock, sch.dt)
        if err < nsolver.min_err
            nsolver.converged = true
            break
        end
        # Check if reaches max_newton_iter
        if nsolver.newton_iter > nsolver.max_newton_iter
            nsolver.converged = false
            break
        end
        newton_step(sim)
    end
    update_dt(sch, fluid, nsolver.converged)
    push!(nsolver.num_iter, nsolver.newton_iter)
    if nsolver.converged
        save_fluid_results(reservoir, sch.t_current)
        save_facility_results(facility, sch.t_current)
        update_fluid_tn(fluid)
    else
        reset_primary_variable(fluid)
    end
    update_phases(fluid, grid.connlist)
    compute_residual(fluid, grid, rock, facility, sch.dt)
    println("NumNewton: $(nsolver.newton_iter)\n")
    nsolver.newton_iter = 0
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

function save_facility_results(facility::Dict{String, AbstractFacility}, t::Float64)
    for w in values(facility)
        save_result(w, t)
    end
end

include("api.jl")
include("owfluid_solver.jl")


end
