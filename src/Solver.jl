module Solver

using SparseArrays:sparse, SparseMatrixCSC

using ..AutoDiff: param, grad, data
using ..Grid: AbstractGrid
using ..State: OWState, get_var_order, M

abstract type NonlinearSolver end

# Newton Raphson Solver
mutable struct NRSolver <: NonlinearSolver
    max_newton_iter::Int
    min_err::Float64
    num_iter::Vector{Int}
    converged::Bool
    δx::Vector{Float64}
    residual::Vector{Float64}
    jac::SparseMatrixCSC{Float64,Int}
    NRSolver() = new(10, 1.0e-6, Vector{Int}(), false)
end

function compute_residual_error(state::OWState, grid::AbstractGrid, dt::Float64)
    a = dt .* M ./ (grid.v .* grid.ϕ)
    rw_err = a .* data(state.bw) .* data(state.rw)
    ro_err = a .* data(state.bo) .* data(state.ro)
    return max(maximum(abs.(rw_err)), maximum((abs.(ro_err))))
end


function assemble_residual(state::OWState)::Vector{Float64}
    numcell = state.numcell
    residual = zeros(Float64, 2*numcell)
    for i = 1:numcell
        residual[2*i-1] = state.rw[i].val
        residual[2*i] = state.ro[i].val
    end
    return residual
end

function assemble_jacobian(state::OWState)::SparseMatrixCSC{Float64,Int}
    I, J, V = Int[], Int[], Float64[]
    rw = state.rw
    ro = state.ro
    numvar = state.numvar
    for i =1:state.numcell
        (ind, g) = grad(rw[i])
        append!(I, (2*i-1) * ones(Int, length(ind) * 2))
        for i = 1:numvar
            append!(J, numvar*ind .+ (i-numvar))
        end
        append!(V, [x for x in transpose(g).data])

        (ind, g) = grad(ro[i])
        append!(I, (2*i) * ones(Int, length(ind) * 2))
        for i = 1:numvar
            append!(J, numvar*ind .+ (i-numvar))
        end
        append!(V, [x for x in transpose(g).data])
    end
    return sparse(I,J,V)
end

function update_solution(state::OWState, δx::Vector{Float64})::Nothing
    nc, nv = state.numcell, state.numvar
    var_order = get_var_order(state)
    ivp, ivsw = var_order["p"], var_order["sw"]
    for i=1:nc
        state.p[i].val -= δx[nv*(i-1) + ivp]
        state.sw[i].val -= δx[nv*(i-1) + ivsw]
    end
    # Bounds check
    for i =1:nc
        if state.sw[i].val < 0.0
            state.sw[i].val = 0.0
        elseif state.sw[i].val > 1.0
            state.sw[i].val = 1.0
        end
    end
end

end
