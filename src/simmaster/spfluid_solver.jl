function compute_residual(
    fluid::SPFluid,
    grid::AbstractGrid,
    rock::AbstractRock,
    wells::Dict{String, AbstractFacility},
    dt::Float64
)
    phase, component = fluid.phases[1], fluid.components[1]
    f = phase.f
    ap, rp = component.a, component.r

    v, ϕ = grid.v, rock.ϕ

    connlist = grid.connlist
    l, r = connlist.l, connlist.r
    # Accumulation term
    @. ap = v * ϕ * (phase.s / phase.b - phase.sn / phase.bn) / dt / M
    @. rp = -value(ap)

    # Flux term
    @views @. rp[l] -= value(f)
    @views @. rp[r] += value(f)

    # Well Term
    for well in values(wells)
        compute_ql(well, fluid)
        # Add well rate to residual
        ind = well.ind
        @. rp[ind] -= value(well.ql)
        @. rp[ind] -= value(well.ql)
    end

    return fluid
end

function assemble_residual(nsolver::NonlinearSolver, fluid::SPFluid)
    @. nsolver.residual = fluid.components[1].r
end

function update_solution(nsolver::NonlinearSolver, fluid::SPFluid)
    phase = fluid.phases[1]
    assembler = nsolver.assembler
    @. assembler.p = value(phase.p) - nsolver.δx
    update_primary_variable(fluid, assembler.p)
    return fluid
end

function compute_residual_error(fluid::SPFluid, grid::AbstractGrid, rock::AbstractRock, dt::Float64)
    a = dt .* M ./ (grid.v .* rock.ϕ)
    p = fluid.phases[1]
    c = fluid.components[1]
    return maximum(abs.(a .* value(p.b) .* c.r))
end


struct SPAssembler <: Assembler
    diag_idx::NamedTuple
    ll_idx::NamedTuple
    rr_idx::NamedTuple
    lr_idx::NamedTuple
    rl_idx::NamedTuple
    p::Vector{Float64}
end

function SPAssembler(
    rowval::Vector{Int},
    colptr::Vector{Int},
    l_list::Vector{Int},
    r_list::Vector{Int},
    nc::Int,
    nconn::Int,
)
    diag_idx = (∂r∂p = Vector{Int64}(undef, nc),)
    for i = 1:nc
        r0 = colptr[i]
        @views rows = rowval[r0:colptr[i+1]-1]
        diag_idx.∂r∂p[i] = r0 + searchsortedfirst(rows, i) - 1
    end
    # Row l, Col l
    ll_idx = (∂f∂p = Vector{Int64}(undef, nconn),)
    # Row r, Col r
    rr_idx = (∂f∂p = Vector{Int64}(undef, nconn),)
    # Row l, Col r, Lower
    lr_idx = (∂f∂p = Vector{Int64}(undef, nconn),)
    # Row r, Col l
    rl_idx = (∂f∂p = Vector{Int64}(undef, nconn),)

    for i = 1:nconn
        l, r = l_list[i], r_list[i]
        # Column l
        r0 = colptr[l]
        @views rows = rowval[r0:colptr[l+1] - 1]
        ll_idx.∂f∂p[i] = r0 + searchsortedfirst(rows, l) - 1
        rl_idx.∂f∂p[i] = r0 + searchsortedfirst(rows, r) - 1
        # Column 2r-1
        r0 = colptr[r]
        @views rows = rowval[r0:colptr[r+1] - 1]
        rr_idx.∂f∂p[i] = r0 + searchsortedfirst(rows, r) - 1
        lr_idx.∂f∂p[i] = r0 + searchsortedfirst(rows, l) - 1
    end
    p = Vector{Float64}(undef, nc)
    return SPAssembler(diag_idx, ll_idx, rr_idx, lr_idx, rl_idx, p)
end


function init_nsolver(nsolver::NonlinearSolver, grid::AbstractGrid, fluid::SPFluid)
    nc = grid.nc
    nsolver.residual = zeros(nc)
    nsolver.δx = zeros(nc)
    ### Init Jacobian
    neighbors, nn = grid.neighbors, grid.num_neighbors
    #
    colptr = Vector{Int}(undef, nc+1)
    colptr[1] = 1
    colptr[2:end] = nn
    cumsum!(colptr, colptr)
    #
    nitem = sum(nn)
    rowval = Vector{Int}(undef, nitem)
    for i in eachindex(nn)
        n, v= nn[i], neighbors[i]
        ofs = colptr[i]
        @. rowval[ofs:ofs+n-1] = v
    end
    nsolver.jac = SparseMatrixCSC{Float64, Int}(nc, nc, colptr, rowval, zeros(nitem))
    # Init Jacobian Assembler
    l_list, r_list, nconn = grid.connlist.l, grid.connlist.r, grid.connlist.nconn
    nsolver.assembler = SPAssembler(rowval, colptr, l_list, r_list, nc, nconn)
end


function assemble_jacobian(nsolver::NonlinearSolver, fluid::SPFluid, wells::Dict{String, AbstractFacility},)
    nzval, assembler = nsolver.jac.nzval, nsolver.assembler
    fill!(nzval, 0.0)
    # Accumulation term
    a = fluid.components[1].a
    diag_idx = assembler.diag_idx
    @. nzval[diag_idx.∂r∂p] = -grad(a, 1, 1)
    # Flux term
    f = fluid.phases[1].f
    ll_idx = assembler.ll_idx
    @views @. nzval[ll_idx.∂f∂p] -= grad(f, 1, 1)
    # Row r, Col l
    rl_idx = assembler.rl_idx
    @views @. nzval[rl_idx.∂f∂p] = grad(f, 1, 1)

    rr_idx = assembler.rr_idx
    @views @. nzval[rr_idx.∂f∂p] += grad(f, 2, 1)

    lr_idx = assembler.lr_idx
    @views @. nzval[lr_idx.∂f∂p] = -grad(f, 2, 1)

    for well in values(wells)
        # Add well rate to residual
        ind = well.ind
        @. nzval[diag_idx.∂r∂p[ind]] -= grad(well.ql, 1, 1)
    end

end
