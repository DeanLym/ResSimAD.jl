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
)
    δx .= jac \ residual
end

lsolver_info(::Julia_BackSlash_Solver) = "Julia backslash"

## BICGSTAB solver with ILU preconditioner
struct BICGSTAB_ILU_Solver <: AbstractLinearSolver
    l::Int
    τ::Float64
    iterations::Vector{Int64}
    BICGSTAB_ILU_Solver(; l=2, τ=0.1) = new(l, τ, Int[])
end

lsolver_info(x::BICGSTAB_ILU_Solver) = "BiCGStab($(x.l)) ILU"

function solve(
    solver::BICGSTAB_ILU_Solver,
    δx::Vector{Float64},
    jac::SparseMatrixCSC{Float64,Int},
    residual::Vector{Float64},
)
    prec = ilu(jac, τ=solver.τ)
    _, log = bicgstabl!(δx, jac, residual, solver.l, Pl=prec, log=true, tol=0.1)
    push!(solver.iterations, log.iters)
end

## GMRES Solver with ILU preconditioner
struct GMRES_ILU_Solver <: AbstractLinearSolver
    τ::Float64
    iterations::Vector{Int64}
    GMRES_ILU_Solver(; τ=0.1) = new(τ, Int[])
end

lsolver_info(::GMRES_ILU_Solver) = "GMRES ILU"

function solve(
    solver::GMRES_ILU_Solver,
    δx::Vector{Float64},
    jac::SparseMatrixCSC{Float64,Int},
    residual::Vector{Float64},
)
    prec = ilu(jac, τ=solver.τ)
    _, log = gmres!(δx, jac, residual, Pl=prec, log=true, tol=0.1)
    push!(solver.iterations, log.iters)
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
    _, log = gmres!(δx, jac, residual, Pl=solver.cpr_prec, log=true, tol=0.1)
    push!(solver.iterations, log.iters)
end

## BICGSTAB Solver with CPR preconditioner
struct BICGSTAB_CPR_Solver <: AbstractLinearSolver
    l::Int
    cpr_prec::CPRPreconditioner
    iterations::Vector{Int64}
end

function BICGSTAB_CPR_Solver(
    nc::Int64,
    neighbors::Vector{Vector{Int}};
    l=2,
    τ=0.5,
)
    cpr_prec = CPRPreconditioner(nc, neighbors, τ)
    return BICGSTAB_CPR_Solver(l, cpr_prec, Int[])
end

lsolver_info(::BICGSTAB_CPR_Solver) = "BICGSTAB CPR"

function solve(
    solver::BICGSTAB_CPR_Solver,
    δx::Vector{Float64},
    jac::SparseMatrixCSC{Float64,Int},
    residual::Vector{Float64},
)::Vector{Float64}
    setup_cpr_preconditioner(solver.cpr_prec, jac)
    _, log = bicgstabl!(δx, jac, residual, solver.l, Pl=solver.cpr_prec, log=true, tol=0.1)
    push!(solver.iterations, log.iters)
end

include("dune_istl_solver.jl")

end
