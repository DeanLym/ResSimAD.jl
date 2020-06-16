module LinearSolver

using SparseArrays:SparseMatrixCSC

using IterativeSolvers, IncompleteLU

abstract type AbstractLinearSolver end

struct Julia_BackSlash_Solver <: AbstractLinearSolver end

struct GMRES_ILU0_Solver <: AbstractLinearSolver
    τ::Float64
    GMRES_ILU0_Solver() = new(0.1)
end

function solve(
    solver::Julia_BackSlash_Solver,
    jac::SparseMatrixCSC{Float64,Int},
    residual::Vector{Float64},
)::Vector{Float64}
    return jac \ residual
end

function solve(
    solver::GMRES_ILU0_Solver,
    jac::SparseMatrixCSC{Float64,Int},
    residual::Vector{Float64},
)::Vector{Float64}
    return gmres(jac, residual, Pl=ilu(jac, τ=solver.τ))
end

end
