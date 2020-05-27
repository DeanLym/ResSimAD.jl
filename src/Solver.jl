module Solver

using SparseArrays:sparse, SparseMatrixCSC

using ..AutoDiff: param, grad
using ..State: OWState, get_var_order

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
    numvar = state.numvar
    I, J, V = Int[], Int[], Float64[]
    for (i, r) in enumerate(state.rw)
        (ind, g) = grad(r)
        for (k, j) in enumerate(ind)
            for iv = 1:numvar
                additem(I,J,V,2*i-1, numvar*(j-1)+iv, g[iv, k])
            end
        end
    end
    for (i, r) in enumerate(state.ro)
        (ind, g) = grad(r)
        for (k, j) in enumerate(ind)
            for iv = 1:numvar
                additem(I,J,V,2*i, numvar*(j-1)+iv, g[iv, k])
            end
        end
    end
    return sparse(I, J, V)
end

function update_solution(state::OWState, δx::Vector{Float64})::Nothing
    nc, nv = state.numcell, state.numvar
    var_order = get_var_order(state)
    ivp, ivsw = var_order["p"], var_order["sw"]
    for i=1:state.numcell
        state.p[i].val -= δx[nv*(i-1) + ivp]
        state.sw[i].val -= δx[nv*(i-1) + ivsw]
    end
end

function additem(
    I::Vector{Int},
    J::Vector{Int},
    V::Vector{Float64},
    i::Int,
    j::Int,
    v::Float64,
)::Nothing
    push!(I, i)
    push!(J, j)
    push!(V, v)
    return nothing
end



end
