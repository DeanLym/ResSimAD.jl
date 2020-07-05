import LinearAlgebra:ldiv!
import Base.\
using AlgebraicMultigrid:ruge_stuben, aspreconditioner, smoothed_aggregation
using AlgebraicMultigrid:Preconditioner

mutable struct CPRPreconditioner
    nc::Int64
    neighbors::Vector{Vector{Int64}}
    𝐑p::Vector{Float64}
    δp::Vector{Float64}
    δx::Vector{Float64}
    𝐉po𝐉oo⁻::Vector{Float64}
    τ::Float64  #ILU fill-in level

    𝐉::SparseMatrixCSC{Float64, Int64}
    amg_prec::Preconditioner
    ilu_prec::IncompleteLU.ILUFactorization{Float64}

    function CPRPreconditioner(nc::Int64, neighbors::Vector{Vector{Int64}}, τ::Float64)
        cpr_prec = new()
        cpr_prec.nc = nc
        cpr_prec.neighbors = neighbors
        cpr_prec.τ = τ
        cpr_prec.𝐑p, cpr_prec.δp = zeros(nc), zeros(nc)
        cpr_prec.δx = zeros(2 * nc)
        cpr_prec.𝐉po𝐉oo⁻ = zeros(nc)
        return cpr_prec
    end
end

function setup_cpr_preconditioner(cpr_prec::CPRPreconditioner, 𝐉::SparseMatrixCSC{Float64, Int64})
    nc, 𝐉po𝐉oo⁻, neighbors = cpr_prec.nc, cpr_prec.𝐉po𝐉oo⁻, cpr_prec.neighbors
    # ILU Preconditioner
    println("Setup ILU")
    @time cpr_prec.ilu_prec = ilu(𝐉, τ=cpr_prec.τ)
    # Construct AMG Preconditioner for pressure equation
    cpr_prec.𝐉 = 𝐉
    # Extract pressure equation (true-IMPES)
    for i = 1:nc
        𝐉po𝐉oo⁻[i] = sum(𝐉[1:2:end, 2*i]) / sum(𝐉[2:2:end, 2*i])
    end
    I, J, V = Int[], Int[], Float64[]
    println("Construct Jp")
    @time for i = 1:nc
        for j in neighbors[i]
            push!(I, i)
            push!(J, j)
            push!(V, 𝐉[2*i-1, 2*j-1] - 𝐉po𝐉oo⁻[i] * 𝐉[2*i, 2*j-1])
        end
    end
    @time 𝐉p = sparse(I, J, V)
    println("AMG Setup")
    @time cpr_prec.amg_prec = aspreconditioner(ruge_stuben(𝐉p))
    # @time cpr_prec.amg_prec = aspreconditioner(smoothed_aggregation(𝐉p))
end


function ldiv!(x::AbstractVector, cpr_prec::CPRPreconditioner, 𝐑::AbstractVector)
    nc, 𝐉po𝐉oo⁻, δx, δp = cpr_prec.nc, cpr_prec.𝐉po𝐉oo⁻, cpr_prec.δx, cpr_prec.δp
    𝐉, 𝐑p = cpr_prec.𝐉 ,cpr_prec.𝐑p
    for i = 1:nc
        𝐑p[i] = 𝐑[2*i-1] - 𝐉po𝐉oo⁻[i] * 𝐑[2*i]
    end
    # println("AMG")
    ldiv!(δp, cpr_prec.amg_prec, 𝐑p) # One AMG iteration for pressure equation
    for i = 1:nc
        δx[2*i-1] = δp[i]
    end
    # println("Jac")
    𝐑2 = 𝐑 .- 𝐉 * δx
    # println("ILU")
    ldiv!(x, cpr_prec.ilu_prec, 𝐑2)
    for i = 1:nc
        x[2*i-1] += δp[i]
    end
    return x
end

function \(cpr_prec::CPRPreconditioner, 𝐑::AbstractVector)
    ldiv!(similar(𝐑), cpr_prec, 𝐑)
end

ldiv!(cpr_prec::CPRPreconditioner, 𝐑::AbstractVector) = copyto!(𝐑, cpr_prec \ 𝐑)
