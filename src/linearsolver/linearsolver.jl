module LinearSolver

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

## GMRES Solver with ILU0 preconditioner
struct GMRES_ILU0_Solver <: AbstractLinearSolver
    τ::Float64
    iterations::Vector{Int64}
    GMRES_ILU0_Solver(; τ=0.1) = new(τ, Int[])
end

function solve(
    solver::GMRES_ILU0_Solver,
    δx::Vector{Float64},
    jac::SparseMatrixCSC{Float64,Int},
    residual::Vector{Float64},
)::Vector{Float64}
    # println("Setup ILU")
    prec = ilu(jac, τ=solver.τ)
    # println("GMRES")
    _, log = gmres!(δx, jac, residual, Pl=prec, log=true)
    # println("Number iteration:", log.iters)
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
    τ=0.1,
)
    cpr_prec = CPRPreconditioner(nc, neighbors, τ)
    return GMRES_CPR_Solver(cpr_prec, Int[])
end

function solve(
    solver::GMRES_CPR_Solver,
    δx::Vector{Float64},
    jac::SparseMatrixCSC{Float64,Int},
    residual::Vector{Float64},
)::Vector{Float64}
    # println("Setup CPR Preconditioner")
    setup_cpr_preconditioner(solver.cpr_prec, jac)
    println("GMRES")
    @time log = gmres!(δx, residual, jac, Pl=solver.cpr_prec, log=true)
    println("Number iteration:", log.iters)
    push!(solver.iterations, log.iters)
    return residual
end

end
