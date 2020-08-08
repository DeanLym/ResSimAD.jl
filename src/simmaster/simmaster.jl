module SimMaster

using LinearAlgebra:norm
using SparseArrays:sparse, SparseMatrixCSC

#! format: off
using ..Global: M, α, β, g_, gc
using ..AutoDiff: advector, value, grad, index
using ..InputParse: parse_input_grid, parse_input_rock, parse_input_fluid

using ..Grid: CartGrid, AbstractGrid, ConnList, construct_connlist, set_cell_depth,
                get_grid_index, set_cell_size, sort_conn
using ..Rock: AbstractRock, StandardRock, set_perm, set_poro
using ..Fluid: AbstractFluid, OWFluid, SPFluid, PVT, PVTC, SWOFTable, SWOFCorey,
                set_fluid_tn, update_phases, update_primary_variable,
                update_fluid_tn, reset_primary_variable
using ..Reservoir: AbstractReservoir, StandardReservoir, save_fluid_results

using ..Facility: AbstractFacility, StandardWell, WellType, Limit, PRODUCER, INJECTOR,
        get_ctrl_mode, compute_wi, compute_qo, compute_qw, compute_ql, compute_bhp,
        save_result

using ..Schedule: Scheduler, update_dt, set_dt, reset_dt, insert_time_step, set_time_step

using ..LinearSolver: AbstractLinearSolver, Julia_BackSlash_Solver,
        GMRES_ILU0_Solver, GMRES_CPR_Solver, solve
#! format: on

abstract type NonlinearSolver end
abstract type Assembler end

@enum Verbose SILENT BRIEF DEBUG ALL

"""
    Sim

Main Sim object
"""
struct Sim
    reservoir::AbstractReservoir
    facility::Dict{String, AbstractFacility}
    scheduler::Scheduler
    nsolver::NonlinearSolver
    lsolver::AbstractLinearSolver
end

mutable struct NRSolver <: NonlinearSolver
    max_newton_iter::Int
    newton_iter::Int
    min_err::Float64
    num_iter::Vector{Int}
    converged::Bool
    recompute_residual::Bool
    δx::Vector{Float64}
    residual::Vector{Float64}
    jac::SparseMatrixCSC{Float64,Int}
    assembler::Assembler
    NRSolver() = new(10, 0, 1.0e-6, Vector{Int}(), false, false)
end

"""
    newton_step(sim::Sim; verbose=BRIEF)

Simulate for one newton step.

# Examples

```jldoctest
julia> using ResSimAD: get_model, newton_step, get_residual_error, SILENT

julia> using Printf: @printf

julia> sim, options = get_model("example1");

julia> @printf("%.3e",get_residual_error(sim))
2.210e-01
julia> newton_step(sim)

julia> @printf("%.3e",get_residual_error(sim))
1.927e-02
julia> newton_step(sim)

julia> @printf("%.3e",get_residual_error(sim))
9.871e-05
```

"""
function newton_step(sim::Sim; verbose=BRIEF)::Nothing
    reservoir = sim.reservoir
    grid, fluid, rock = reservoir.grid, reservoir.fluid, reservoir.rock
    facility, nsolver, lsolver, sch = sim.facility, sim.nsolver, sim.lsolver, sim.scheduler
    # Assemble residual
    assemble_residual(nsolver, fluid)
    # Assemble jacobian
    # println("\nAssemble Jacobian")
    assemble_jacobian(nsolver, fluid, facility)
    # Solve equation
    # println("Solve Equation")
    solve(lsolver, nsolver.δx, nsolver.jac, nsolver.residual)
    # println("Solved Equation")
    # Update solution
    update_solution(nsolver, fluid)
    # Update dynamic states and then compute new residual
    # println("Update Phases")
    update_phases(fluid, grid.connlist)
    # println("Compute Residual")
    compute_residual(fluid, grid, rock, facility, sch.dt)
    nsolver.newton_iter += 1
    if verbose >= DEBUG println("Newton Step ", nsolver.newton_iter) end
    return nothing
end

"""
    time_step(sim::Sim; verbose=BRIEF)

Simluate for one time step. Time step length is `sim.scheduler.dt`.

# Examples

```jldoctest
julia> using ResSimAD: get_model, time_step, SILENT

julia> sim, options = get_model("example1");

julia> sim.scheduler.t_current
0.0

julia> sim.scheduler.dt
0.1

julia> time_step(sim, verbose=SILENT);

julia> sim.scheduler.t_current
0.1
```

"""
function time_step(sim::Sim; verbose=BRIEF)::Nothing
    reservoir = sim.reservoir
    grid, fluid, rock = reservoir.grid, reservoir.fluid, reservoir.rock
    facility, nsolver, lsolver, sch = sim.facility, sim.nsolver, sim.lsolver, sim.scheduler
    if verbose >= BRIEF println("Day ", sch.t_next) end
    while true
        if nsolver.recompute_residual
            compute_residual(fluid, grid, rock, facility, sch.dt)
            nsolver.recompute_residual = false
        end
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
        newton_step(sim; verbose=verbose)
    end
    update_dt(sch, fluid, nsolver.converged)
    push!(nsolver.num_iter, nsolver.newton_iter)
    if nsolver.converged
        save_fluid_results(reservoir, sch.t_current)
        save_facility_results(facility, sch.t_current)
        update_fluid_tn(fluid)
    else
        reset_primary_variable(fluid)
        if verbose >= BRIEF println("[WARNING] Converge failed, cut time step.") end
    end
    update_phases(fluid, grid.connlist)
    compute_residual(fluid, grid, rock, facility, sch.dt)
    if verbose >= BRIEF println("  NumNewton:", nsolver.newton_iter) end
    nsolver.newton_iter = 0
    return nothing
end

"""
    step_to(sim::Sim, t::Float64; verbose=BRIEF)

Simluate until day `t`.

# Examples
```jldoctest
julia> using ResSimAD: get_model, step_to, SILENT

julia> sim, options = get_model("example1");

julia> sim.scheduler.t_current
0.0

julia> step_to(sim, 100.0, verbose=SILENT);

julia> sim.scheduler.t_current
100.0
```

"""
function step_to(sim::Sim, t::Float64; verbose=BRIEF)::Nothing
    sch = sim.scheduler
    insert_time_step(sch, t)
    while sch.t_current < t
        time_step(sim; verbose=verbose)
    end
end

"""
    runsim(sim::Sim; verbose=BRIEF)

Run simulation from day `sim.scheduler.t_current` to day `sim.scheduler.time_step[end]`

# Examples
```jldoctest
julia> using ResSimAD: get_model, runsim, SILENT

julia> sim, options = get_model("example1");

julia> sim.scheduler.t_current
0.0

julia> sim.scheduler.time_step[end]
1825.0

julia> runsim(sim, verbose=SILENT);

julia> sim.scheduler.t_current
1825.0
```

"""
function runsim(sim::Sim; verbose=BRIEF)::Nothing
    sch = sim.scheduler
    while sch.t_current < sch.time_step[end]
        time_step(sim; verbose=verbose)
    end
    return nothing
end


function get_residual_error(sim::Sim)
    reservoir = sim.reservoir
    grid, fluid, rock = reservoir.grid, reservoir.fluid, reservoir.rock
    sch = sim.scheduler
    compute_residual_error(fluid, grid, rock, sch.dt)
end

function save_facility_results(facility::Dict{String, AbstractFacility}, t::Float64)
    for w in values(facility)
        save_result(w, t)
    end
end

include("api.jl")
include("owfluid_solver.jl")
include("spfluid_solver.jl")

end
