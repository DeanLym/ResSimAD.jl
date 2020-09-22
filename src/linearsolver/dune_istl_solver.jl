using Libdl
using ..Grid: AbstractGrid
using duneistl_bicgilu_julia_jll

const dune_istl_lib_path = duneistl_bicgilu_julia_jll.libduneistl_bicgilu

## BICGSTAB solver with ILU preconditioner from Dune-ISTL
struct BICGSTAB_ILU_DUNE_ISTL_Solver <: AbstractLinearSolver
    lib_ptr::Ptr{nothing}
    iterations::Vector{Int64}
    BI::Vector{Int64}
    BJ::Vector{Int64}
    row_size::Vector{Int64}
    function BICGSTAB_ILU_DUNE_ISTL_Solver(grid::AbstractGrid)
        lib_ptr = dlopen(dune_istl_lib_path)
        BI = Int[]
        BJ = Int[]
        row_size = Int[]
        for (i, nn) in enumerate(grid.neighbors)
            for j in nn
                push!(BI, i)
                push!(BJ, j)
            end
        end
        for n in grid.num_neighbors
            push!(row_size, n)
        end
        return new(lib_ptr, Int[], BI, BJ, row_size)
    end
end

lsolver_info(::BICGSTAB_ILU_DUNE_ISTL_Solver) = "BiCGStab ILU (DUNE-ISTL)"


function solve(
    solver::BICGSTAB_ILU_DUNE_ISTL_Solver,
    δx::Vector{Float64},
    jac::SparseMatrixCSC{Float64,Int},
    residual::Vector{Float64},
)
    V1 = Vector{Float64}()
    V2 = Vector{Float64}()
    V3 = Vector{Float64}()
    V4 = Vector{Float64}()

    BI, BJ = solver.BI, solver.BJ

    for (i, j) in zip(BI, BJ)
        push!(V1, jac[2i-1, 2j-1])
        push!(V2, jac[2i-1, 2j])
        push!(V3, jac[2i, 2j-1])
        push!(V4, jac[2i, 2j])
    end

    n = size(jac)[1] ÷ 2 
    nnz = length(V1)

    δx1 = δx[1:2:end]
    δx2 = δx[2:2:end]

    tol = 1.0e-2
    max_iter = 100

    ccall(dlsym(solver.lib_ptr, "istl_solve_block"), Cint, 
    (Cint, Cint, Ref{Cint}, Ref{Cint}, Ref{Cint}, 
    Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},
    Ref{Cdouble},Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Cdouble, Cint), 
    Int32(n), Int32(nnz), Int32.(solver.row_size), Int32.(BI .- 1), Int32.(BJ .- 1), 
    V1, V2, V3, V4, residual[1:2:end], residual[2:2:end], δx1, δx2, tol, max_iter) 

    @. δx[1:2:end] = δx1
    @. δx[2:2:end] = δx2

    push!(solver.iterations, 0)

end
