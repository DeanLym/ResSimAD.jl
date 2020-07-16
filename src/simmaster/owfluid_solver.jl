function compute_residual(
    fluid::OWFluid,
    grid::AbstractGrid,
    rock::AbstractRock,
    wells::Dict{String, AbstractFacility},
    dt::Float64
)
    oc, wc = fluid.components.o, fluid.components.w
    o, w = fluid.phases.o, fluid.phases.w
    v, ϕ = grid.v, rock.ϕ
    ao, aw = oc.a, wc.a
    fo, fw = o.f, w.f
    ro, rw = oc.r, wc.r
    connlist = grid.connlist
    l, r = connlist.l, connlist.r
    # Accumulation term
    @. ao = v * ϕ * (o.s / o.b - o.sn / o.bn) / dt / M
    @. aw = v * ϕ * (w.s / w.b - w.sn / w.bn) / dt / M

    @. ro = -value(ao)
    @. rw = -value(aw)

    # Flux term
    @views @. ro[l] -= value(fo)
    @views @. ro[r] += value(fo)
    @views @. rw[l] -= value(fw)
    @views @. rw[r] += value(fw)

    # Well Term
    for well in values(wells)
        compute_qo(well, fluid)
        compute_qw(well, fluid)
        # Add well rate to residual
        ind = well.ind
        @. ro[ind] -= value(well.qo)
        @. rw[ind] -= value(well.qw)
    end

    return fluid
end

function assemble_residual(nsolver::NonlinearSolver, fluid::OWFluid)
    @. nsolver.residual[1:2:end] = fluid.components.w.r
    @. nsolver.residual[2:2:end] = fluid.components.o.r
end

function update_solution(nsolver::NonlinearSolver, fluid::OWFluid)
    o, w = fluid.phases.o, fluid.phases.w
    assembler = nsolver.assembler
    @. assembler.po = value(o.p) - nsolver.δx[1:2:end]
    @. assembler.sw = value(w.s) - nsolver.δx[2:2:end]
    # Bounds check
    @. assembler.sw = max(assembler.sw, 0.0)
    @. assembler.sw = min(assembler.sw, 1.0)
    update_primary_variable(fluid, assembler.po, assembler.sw)
    return fluid
end


function compute_residual_error(fluid::OWFluid, grid::AbstractGrid, rock::AbstractRock, dt::Float64)
    a = dt .* M ./ (grid.v .* rock.ϕ)
    op, wp = fluid.phases.o, fluid.phases.w
    oc, wc = fluid.components.o, fluid.components.w
    rw_err = maximum(abs.(a .* value(wp.b) .* wc.r))
    ro_err = maximum(abs.(a .* value(op.b) .* oc.r))
    return max(rw_err, ro_err)
end

struct OWAssembler <: Assembler
    diag_idx::NamedTuple
    ll_idx::NamedTuple
    rr_idx::NamedTuple
    lr_idx::NamedTuple
    rl_idx::NamedTuple
    po::Vector{Float64}
    sw::Vector{Float64}
end


function OWAssembler(
    rowval::Vector{Int},
    colptr::Vector{Int},
    l_list::Vector{Int},
    r_list::Vector{Int},
    nc::Int,
    nconn::Int,
)
    diag_idx = (∂rw∂p = Vector{Int64}(undef, nc),
                ∂rw∂s = Vector{Int64}(undef, nc),
                ∂ro∂p = Vector{Int64}(undef, nc),
                ∂ro∂s = Vector{Int64}(undef, nc),)
    for i = 1:nc
        r0 = colptr[2i-1]
        @views rows = rowval[r0:colptr[2i]-1]
        diag_idx.∂rw∂p[i] = r0 + searchsortedfirst(rows, 2i-1) - 1
        diag_idx.∂ro∂p[i] = r0 + searchsortedfirst(rows, 2i) - 1

        r0 = colptr[2i]
        @views rows = rowval[r0:colptr[2i+1]-1]
        diag_idx.∂rw∂s[i] = r0 + searchsortedfirst(rows, 2i-1) - 1
        diag_idx.∂ro∂s[i] = r0 + searchsortedfirst(rows, 2i) - 1
    end
    # Row l, Col l
    ll_idx = (∂fw∂p = Vector{Int64}(undef, nconn),
              ∂fw∂s = Vector{Int64}(undef, nconn),
              ∂fo∂p = Vector{Int64}(undef, nconn),
              ∂fo∂s = Vector{Int64}(undef, nconn),)
    # Row r, Col r
    rr_idx = (∂fw∂p = Vector{Int64}(undef, nconn),
              ∂fw∂s = Vector{Int64}(undef, nconn),
              ∂fo∂p = Vector{Int64}(undef, nconn),
              ∂fo∂s = Vector{Int64}(undef, nconn),)
    # Row l, Col r, Lower
    lr_idx = (∂fw∂p = Vector{Int64}(undef, nconn),
              ∂fw∂s = Vector{Int64}(undef, nconn),
              ∂fo∂p = Vector{Int64}(undef, nconn),
              ∂fo∂s = Vector{Int64}(undef, nconn),)
    # Row r, Col l
    rl_idx = (∂fw∂p = Vector{Int64}(undef, nconn),
              ∂fw∂s = Vector{Int64}(undef, nconn),
              ∂fo∂p = Vector{Int64}(undef, nconn),
              ∂fo∂s = Vector{Int64}(undef, nconn),)

    for i = 1:nconn
        l, r = l_list[i], r_list[i]
        # Column 2l-1
        r0 = colptr[2l-1]
        @views rows = rowval[r0:colptr[2l] - 1]
        ll_idx.∂fw∂p[i] = r0 + searchsortedfirst(rows, 2l-1) - 1
        ll_idx.∂fo∂p[i] = r0 + searchsortedfirst(rows, 2l) - 1
        rl_idx.∂fw∂p[i] = r0 + searchsortedfirst(rows, 2r-1) - 1
        rl_idx.∂fo∂p[i] = r0 + searchsortedfirst(rows, 2r) - 1
        # Column 2l
        r0 = colptr[2l]
        @views rows = rowval[r0:colptr[2l+1] - 1]
        ll_idx.∂fw∂s[i] = r0 + searchsortedfirst(rows, 2l-1) - 1
        ll_idx.∂fo∂s[i] = r0 + searchsortedfirst(rows, 2l) - 1
        rl_idx.∂fw∂s[i] = r0 + searchsortedfirst(rows, 2r-1) - 1
        rl_idx.∂fo∂s[i] = r0 + searchsortedfirst(rows, 2r) - 1
        # Column 2r-1
        r0 = colptr[2r-1]
        @views rows = rowval[r0:colptr[2r] - 1]
        rr_idx.∂fw∂p[i] = r0 + searchsortedfirst(rows, 2r-1) - 1
        rr_idx.∂fo∂p[i] = r0 + searchsortedfirst(rows, 2r) - 1
        lr_idx.∂fw∂p[i] = r0 + searchsortedfirst(rows, 2l-1) - 1
        lr_idx.∂fo∂p[i] = r0 + searchsortedfirst(rows, 2l) - 1
        # Column 2r
        r0 = colptr[2r]
        @views rows = rowval[r0:colptr[2r+1] - 1]
        rr_idx.∂fw∂s[i] = r0 + searchsortedfirst(rows, 2r-1) - 1
        rr_idx.∂fo∂s[i] = r0 + searchsortedfirst(rows, 2r) - 1
        lr_idx.∂fw∂s[i] = r0 + searchsortedfirst(rows, 2l-1) - 1
        lr_idx.∂fo∂s[i] = r0 + searchsortedfirst(rows, 2l) - 1
    end
    po = Vector{Float64}(undef, nc)
    sw = Vector{Float64}(undef, nc)
    return OWAssembler(diag_idx, ll_idx, rr_idx, lr_idx, rl_idx, po, sw)
end


function init_nsolver(nsolver::NonlinearSolver, grid::AbstractGrid, fluid::OWFluid)
    nc = grid.nc
    nsolver.residual = zeros(2*nc)
    nsolver.δx = zeros(2*nc)
    ### Init Jacobian
    neighbors, nn = grid.neighbors, grid.num_neighbors
    #
    colptr = Vector{Int}(undef, 2nc+1)
    colptr[1] = 1
    colptr[2:2:end] = 2*nn
    colptr[3:2:end] = 2*nn
    cumsum!(colptr, colptr)
    #
    nitem = 4*sum(nn)
    rowval = Vector{Int}(undef, nitem)
    for i in eachindex(nn)
        n, v= nn[i], neighbors[i]
        ofs = colptr[2i-1]
        @. rowval[ofs:2:ofs+2n-1] = 2v-1
        @. rowval[ofs+1:2:ofs+2n] = 2v
        ofs = colptr[2i]
        @. rowval[ofs:2:ofs+2n-1] = 2v-1
        @. rowval[ofs+1:2:ofs+2n] = 2v
    end
    nsolver.jac = SparseMatrixCSC{Float64, Int}(2nc, 2nc, colptr, rowval, zeros(nitem))
    # Init Jacobian Assembler
    l_list, r_list, nconn = grid.connlist.l, grid.connlist.r, grid.connlist.nconn
    nsolver.assembler = OWAssembler(rowval, colptr, l_list, r_list, nc, nconn)
end


function assemble_jacobian(nsolver::NonlinearSolver, fluid::OWFluid, wells::Dict{String, AbstractFacility},)
    nzval, assembler = nsolver.jac.nzval, nsolver.assembler
    # Accumulation term
    aw, ao = fluid.components.w.a, fluid.components.o.a
    diag_idx = assembler.diag_idx
    @. nzval[diag_idx.∂rw∂p] = -grad(aw, 1, 1)
    @. nzval[diag_idx.∂ro∂p] = -grad(ao, 1, 1)
    @. nzval[diag_idx.∂rw∂s] = -grad(aw, 1, 2)
    @. nzval[diag_idx.∂ro∂s] = -grad(ao, 1, 2)
    # Flux term
    fw, fo = fluid.phases.w.f, fluid.phases.o.f
    ll_idx = assembler.ll_idx
    @views @. nzval[ll_idx.∂fw∂p] -= grad(fw, 1, 1)
    @views @. nzval[ll_idx.∂fo∂p] -= grad(fo, 1, 1)
    @views @. nzval[ll_idx.∂fw∂s] -= grad(fw, 1, 2)
    @views @. nzval[ll_idx.∂fo∂s] -= grad(fo, 1, 2)
    # Row r, Col l
    rl_idx = assembler.rl_idx
    @views @. nzval[rl_idx.∂fw∂p] = grad(fw, 1, 1)
    @views @. nzval[rl_idx.∂fo∂p] = grad(fo, 1, 1)
    @views @. nzval[rl_idx.∂fw∂s] = grad(fw, 1, 2)
    @views @. nzval[rl_idx.∂fo∂s] = grad(fo, 1, 2)

    rr_idx = assembler.rr_idx
    @views @. nzval[rr_idx.∂fw∂p] += grad(fw, 2, 1)
    @views @. nzval[rr_idx.∂fo∂p] += grad(fo, 2, 1)
    @views @. nzval[rr_idx.∂fw∂s] += grad(fw, 2, 2)
    @views @. nzval[rr_idx.∂fo∂s] += grad(fo, 2, 2)

    lr_idx = assembler.lr_idx
    @views @. nzval[lr_idx.∂fw∂p] = -grad(fw, 2, 1)
    @views @. nzval[lr_idx.∂fo∂p] = -grad(fo, 2, 1)
    @views @. nzval[lr_idx.∂fw∂s] = -grad(fw, 2, 2)
    @views @. nzval[lr_idx.∂fo∂s] = -grad(fo, 2, 2)

    for well in values(wells)
        # Add well rate to residual
        ind = well.ind
        @. nzval[diag_idx.∂rw∂p[ind]] -= grad(well.qw, 1, 1)
        @. nzval[diag_idx.∂rw∂s[ind]] -= grad(well.qw, 1, 2)
        @. nzval[diag_idx.∂ro∂p[ind]] -= grad(well.qo, 1, 1)
        @. nzval[diag_idx.∂ro∂s[ind]] -= grad(well.qo, 1, 2)
    end
end
