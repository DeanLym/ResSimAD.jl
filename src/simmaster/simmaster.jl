module SimMaster

using CSV
using HDF5
using Memento
using Formatting:format, FormatExpr, fmt, FormatSpec

const LOGGER = getlogger(@__MODULE__)
__init__() = Memento.register(LOGGER)

using LinearAlgebra:norm
using SparseArrays:sparse, SparseMatrixCSC
using Statistics
using DuneIstlSolvers

#! format: off
using ..Global: M, α, β, g_, gc

using ..AutoDiff: advector, value, grad, index

using ..InputParser: parse_input, parse_well_option

using ..Grid: CartGrid, AbstractGrid, ConnList, construct_connlist, set_cell_depth,
                get_grid_index, set_cell_size, sort_conn, grid_info, compute_cell_depth

using ..Rock: AbstractRock, StandardRock, TransRock, set_perm, set_poro, set_trans

using ..Fluid: AbstractFluid, OWFluid, SPFluid, PVT, PVTC, SWOFTable, SWOFCorey,
                set_fluid_tn, update_phases, update_primary_variable,
                update_fluid_tn, reset_primary_variable, fluid_system, equil,
                primary_variables

using ..Reservoir: AbstractReservoir, StandardReservoir, save_fluid_results

using ..Facility: AbstractFacility, StandardWell, WellType, Limit, PRODUCER, INJECTOR,
        get_ctrl_mode, get_limit, compute_wi, compute_well_state, save_result, check_limits

using ..Schedule: Scheduler, update_dt, set_dt, reset_dt, insert_time_step, set_time_step

using ..LinearSolver: AbstractLinearSolver, Julia_BackSlash_Solver, DuneIstlSolverWrapper,
        GMRES_ILU_Solver, GMRES_CPR_Solver, BICGSTAB_ILU_Solver, BICGSTAB_CPR_Solver, lsolve, lsolver_info
#! format: on

import Base.show

abstract type AbstractNonlinearSolver end
abstract type AbstractAssembler end

const log_fmt = FormatExpr("Day = {:>10.3f}, ΔT = {:>8.3f}, NI = {:>4d}, LI ={:>6d}, t = {:>8.3f}secs")

"""
    Sim

Main Sim object
"""
 mutable struct Sim
    reservoir::AbstractReservoir
    facility::Dict{String, AbstractFacility}
    scheduler::Scheduler
    nsolver::AbstractNonlinearSolver
end

function Base.show(io::IO, sim::Sim)
    println(io, "Sim model")
    println(io, "Grid: ", grid_info(sim.reservoir.grid))
    println(io, "Number of cells: ", sim.nc)
    println(io, "Number of connections: ", sim.connlist.nconn)
    println(io, "Fluid system: ", fluid_system(sim.reservoir.fluid))
    println(io, "Wells: ", collect(keys(sim.facility)))
    println(io, "Linear solver: ", lsolver_info(sim.nsolver.lsolver))
end


mutable struct NRSolverDuneIstl <: AbstractNonlinearSolver
    max_newton_iter::Int
    newton_iter::Int
    min_err::Float64
    num_iter::Vector{Int}
    converged::Bool
    recompute_residual::Bool
    lsolver::DuneIstlSolverWrapper
    assembler::AbstractAssembler
    NRSolverDuneIstl() = new(10, 0, 1.0e-6, Vector{Int}(), false, false)
end

function solve_linear_equations(nsolver::NRSolverDuneIstl)
    lsolve(nsolver.lsolver)
end

mutable struct NRSolver <: AbstractNonlinearSolver
    max_newton_iter::Int
    newton_iter::Int
    min_err::Float64
    num_iter::Vector{Int}
    converged::Bool
    recompute_residual::Bool
    δx::Vector{Float64}
    residual::Vector{Float64}
    jac::SparseMatrixCSC{Float64,Int}
    assembler::AbstractAssembler
    lsolver::AbstractLinearSolver
    NRSolver() = new(10, 0, 1.0e-6, Vector{Int}(), false, false)
end

function solve_linear_equations(nsolver::NRSolver)
    lsolve(nsolver.lsolver, nsolver.δx, nsolver.jac, nsolver.residual)
end

"""
    newton_step(sim::Sim)

Simulate for one newton step.

# Examples

```jldoctest
julia> using ResSimAD: get_model, newton_step, get_residual_error

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
function newton_step(sim::Sim)::Nothing
    reservoir = sim.reservoir
    grid, fluid, rock = reservoir.grid, reservoir.fluid, reservoir.rock
    facility, nsolver, sch = sim.facility, sim.nsolver, sim.scheduler
    # Assemble residual
    assemble_residual(nsolver, fluid)
    # Assemble jacobian
    assemble_jacobian(nsolver, fluid, facility)
    # Solve equation
    solve_linear_equations(nsolver)
    # Update solution
    update_solution(nsolver, fluid)
    # Update dynamic states and then compute new residual
    update_phases(fluid, grid.connlist)
    compute_residual(fluid, grid, rock, facility, sch.dt)
    nsolver.newton_iter += 1
    return nothing
end


"""
    time_step(sim::Sim)

Simluate for one time step. Time step length is `sim.scheduler.dt`.

# Examples

```jldoctest
julia> using ResSimAD: get_model, time_step

julia> sim, options = get_model("example1");

julia> sim.scheduler.t_current
0.0

julia> sim.scheduler.dt
0.1

julia> time_step(sim);

julia> sim.scheduler.t_current
0.1
```

"""
function time_step(sim::Sim)::Nothing
    t0 = time()
    reservoir = sim.reservoir
    grid, fluid, rock = reservoir.grid, reservoir.fluid, reservoir.rock
    facility, nsolver, lsolver, sch = sim.facility, sim.nsolver, sim.nsolver.lsolver, sim.scheduler
    dt, t = sch.dt, sch.t_next
    linear_iter = 0
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
        newton_step(sim)
        linear_iter += lsolver.iterations[end]
    end
    push!(nsolver.num_iter, nsolver.newton_iter)

    proceed = true
    if nsolver.converged
        for well in values(facility)
            if !(check_limits(well))
                proceed = false
            end
        end
    else
        proceed = false
    end
    
    if proceed
        update_dt(sch, fluid, nsolver.converged)
        save_fluid_results(reservoir, sch.t_current)
        save_facility_results(facility, sch.t_current)
        update_fluid_tn(fluid)
    else
        if !nsolver.converged
            update_dt(sch, fluid, nsolver.converged)
            warn(LOGGER, "Convergence failed, cutting time step size to $(round(sch.dt, digits=3)) days")
        end
        reset_primary_variable(fluid)
    end

    update_phases(fluid, grid.connlist)
    compute_residual(fluid, grid, rock, facility, sch.dt)
    if proceed
        info(LOGGER, format(log_fmt, t, dt, nsolver.newton_iter, linear_iter, time() - t0))
    end
    nsolver.newton_iter = 0
    return nothing
end

"""
    step_to(sim::Sim, t::Float64)

Simluate until day `t`.

# Examples
```jldoctest
julia> using ResSimAD: get_model, step_to

julia> sim, options = get_model("example1");

julia> sim.scheduler.t_current
0.0

julia> step_to(sim, 100.0);

julia> sim.scheduler.t_current
100.0
```

"""
function step_to(sim::Sim, t::Float64)::Nothing
    sch = sim.scheduler
    insert_time_step(sch, t)
    while sch.t_current < t
        time_step(sim)
    end
end

"""
    runsim(sim::Sim)

Run simulation from day `sim.scheduler.t_current` to day `sim.scheduler.time_step[end]`

# Examples
```jldoctest
julia> using ResSimAD: get_model, runsim

julia> sim, options = get_model("example1");

julia> sim.scheduler.t_current
0.0

julia> sim.scheduler.time_step[end]
1825.0

julia> runsim(sim);

julia> sim.scheduler.t_current
1825.0
```

"""
function runsim(sim::Sim)::Nothing
    sch = sim.scheduler
    t0 = time()
    while sch.t_current < sch.time_step[end]
        time_step(sim)
    end
    dt = round(time() - t0; digits = 3)
    info(LOGGER, "Simulation finished, elapsed time: $(dt) seconds")
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

include("sim_constructor.jl")
include("api_functions.jl")
include("derived_property.jl")

include("owfluid_solver.jl")
include("spfluid_solver.jl")



end
