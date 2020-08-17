module LinearSolver

using Memento

const LOGGER = getlogger(@__MODULE__)
__init__() = Memento.register(LOGGER)

using SparseArrays:SparseMatrixCSC, sparse
using IterativeSolvers, IncompleteLU


abstract type AbstractLinearSolver end

## Julia BackSlash Solver
struct Julia_BackSlash_Solver <: AbstractLinearSolver end

function solve(
    solver::Julia_BackSlash_Solver,
    δx::Vector{Float64},
    jac::SparseMatrixCSC{Float64,Int},
    residual::Vector{Float64},
)::Vector{Float64}
    residual .= jac \ residual
end

lsolver_info(::Julia_BackSlash_Solver) = "Julia backslash"

## GMRES Solver with ILU0 preconditioner
struct GMRES_ILU0_Solver <: AbstractLinearSolver
    τ::Float64
    iterations::Vector{Int64}
    GMRES_ILU0_Solver(; τ=0.1) = new(τ, Int[])
end

lsolver_info(::GMRES_ILU0_Solver) = "GMRES ILU"

function solve(
    solver::GMRES_ILU0_Solver,
    δx::Vector{Float64},
    jac::SparseMatrixCSC{Float64,Int},
    residual::Vector{Float64},
)::Vector{Float64}
    prec = ilu(jac, τ=solver.τ)
    _, log = gmres!(δx, jac, residual, Pl=prec, log=true)
    push!(solver.iterations, log.iters)
    return residual
end

## GMRES Solver with CPR preconditioner

include("cpr_preconditioner.jl")

struct GMRES_CPR_Solver <: AbstractLinearSolver
    cpr_prec::CPRPreconditioner
    iterations::Vector{Int64}
end

function GMRES_CPR_Solver(
    nc::Int64,
    neighbors::Vector{Vector{Int}};
    τ=0.5,
)
    cpr_prec = CPRPreconditioner(nc, neighbors, τ)
    return GMRES_CPR_Solver(cpr_prec, Int[])
end

lsolver_info(::GMRES_CPR_Solver) = "GMRES CPR"

function solve(
    solver::GMRES_CPR_Solver,
    δx::Vector{Float64},
    jac::SparseMatrixCSC{Float64,Int},
    residual::Vector{Float64},
)::Vector{Float64}
    setup_cpr_preconditioner(solver.cpr_prec, jac)
    # _, log = gmres!(δx, jac, residual, Pl=prec, log=true)
    _, log = gmres!(δx, jac, residual, Pl=solver.cpr_prec, log=true, maxiter=200)
    push!(solver.iterations, log.iters)
    return residual
end

end
