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
        compute_well_state(well, fluid)
        # Add well rate to residual
        ind = well.ind
        @. ro[ind] -= value(well.qo)
        @. rw[ind] -= value(well.qw)
    end

    return fluid
end

function assemble_residual(nsolver::NRSolver, fluid::OWFluid)
    @. nsolver.residual[1:2:end] = fluid.components.w.r
    @. nsolver.residual[2:2:end] = fluid.components.o.r
end

function update_solution(nsolver::NRSolver, fluid::OWFluid)
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

function compute_residual_error(fluid::OWFluid, grid::AbstractGrid, rock::AbstractRock, wells::Dict{String, AbstractFacility}, dt::Float64)
    a = dt .* M ./ (grid.v .* rock.ϕ)
    op, wp = fluid.phases.o, fluid.phases.w
    oc, wc = fluid.components.o, fluid.components.w
    rw_err = maximum(abs.(a .* value(wp.b) .* wc.r))
    ro_err = maximum(abs.(a .* value(op.b) .* oc.r))
    w_err = maximum(abs.([w.rw for w in values(wells)]))
    return max(rw_err, ro_err, w_err)
end

struct OWAssembler <: AbstractAssembler
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


function initialize_nsolver(nsolver::NRSolver, grid::AbstractGrid, wells::Dict{String, AbstractFacility}, ::OWFluid)
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


function assemble_jacobian(nsolver::NRSolver, fluid::OWFluid, wells::Dict{String, AbstractFacility},)
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


## DUNE-ISTL Backend
struct OWAssemblerDuneIstl <: AbstractAssembler
    nc::Int
    nconn::Int
    diag_bi::Vector{Int32}
    l_bi::Vector{Int32}
    r_bi::Vector{Int32}
    v_nc::Vector{Float64}
    v_nconn::Vector{Float64}
    po::Vector{Float64}
    sw::Vector{Float64}
end

function OWAssemblerDuneIstl(
    grid::AbstractGrid,
)
    nc, nconn = grid.nc, grid.connlist.nconn
    diag_bi = collect(Int32, 1:nc)
    l_bi, r_bi = Int32.(grid.connlist.l), Int32.(grid.connlist.r)
    v_nc = Vector{Float64}(undef, nc)
    v_nconn = Vector{Float64}(undef, nconn)
    po = Vector{Float64}(undef, nc)
    sw = Vector{Float64}(undef, nc)
    return OWAssemblerDuneIstl(nc, nconn, diag_bi, l_bi, r_bi, v_nc, v_nconn, po, sw)
end


function initialize_nsolver(nsolver::NRSolverDuneIstl, grid::AbstractGrid, wells::Dict{String, AbstractFacility}, ::AbstractFluid)
    indices = Set{Tuple{Int, Int}}()
    row_size = Int[]
    for (i, nn) in enumerate(grid.neighbors)
        for j in nn
            push!(indices, (i, j))
        end
    end
    for n in grid.num_neighbors
        push!(row_size, n)
    end
    # Add well block
    for w in values(wells)
        for i in w.ind
            for j in w.ind
                if !((i,j)∈indices)
                    push!(indices, (i, j))
                    row_size[i] += 1
                end
            end
        end
    end
    BI = [x[1] for x in indices]
    BJ = [x[2] for x in indices]

    nnz = length(BI)

    construct_matrix(nsolver.lsolver.solver, nnz, Int32.(row_size), Int32.(BI), Int32.(BJ))
    nsolver.assembler = OWAssemblerDuneIstl(grid)
end

function assemble_residual(nsolver::NRSolverDuneIstl, fluid::OWFluid, wells::Dict{String, AbstractFacility})
    solver, assembler = nsolver.lsolver.solver, nsolver.assembler
    n = assembler.nc
    reset_rhs(solver)
    add_value_rhs(solver, n, assembler.diag_bi, 1, fluid.components.w.r)
    add_value_rhs(solver, n, assembler.diag_bi, 2, fluid.components.o.r)
    # Add well terms
    for well in values(wells)
        ind = Int32.(well.ind)
        nperf = length(ind)
        ∂rw∂w, ∂ro∂w = -grad(well.qw, 2, 1), -grad(well.qo, 2, 1) # Cw
        ∂w∂w = well.∂w # Dw
        Rw = well.rw
        add_value_rhs(solver, nperf, ind, 1, - ∂rw∂w .* Rw ./ ∂w∂w)
        add_value_rhs(solver, nperf, ind, 2, - ∂ro∂w .* Rw ./ ∂w∂w)
    end
end

function assemble_jacobian(nsolver::NRSolverDuneIstl, fluid::OWFluid, wells::Dict{String, AbstractFacility})
    solver, assembler = nsolver.lsolver.solver, nsolver.assembler
    nc, nconn = assembler.nc, assembler.nconn
    reset_matrix(solver, false)
    # Accumulation term
    aw, ao = fluid.components.w.a, fluid.components.o.a
    diag_bi = assembler.diag_bi
    v_nc, v_nconn = assembler.v_nc, assembler.v_nconn
    # void add_value_matrix(int nn, int* BI, int* BJ, int I, int J, T* value){
    @. v_nc = -grad(aw, 1, 1)
    add_value_matrix(solver, nc, diag_bi, diag_bi, 1, 1, v_nc)
    @. v_nc = -grad(ao, 1, 1)
    add_value_matrix(solver, nc, diag_bi, diag_bi, 2, 1, v_nc)
    @. v_nc = -grad(aw, 1, 2)
    add_value_matrix(solver, nc, diag_bi, diag_bi, 1, 2, v_nc)
    @. v_nc = -grad(ao, 1, 2)
    add_value_matrix(solver, nc, diag_bi, diag_bi, 2, 2, v_nc)
    # # Flux term
    fw, fo = fluid.phases.w.f, fluid.phases.o.f
    l_bi, r_bi = assembler.l_bi, assembler.r_bi
    @. v_nconn = -grad(fw, 1, 1)
    add_value_matrix(solver, nconn, l_bi, l_bi, 1, 1, v_nconn)
    @. v_nconn = -grad(fo, 1, 1)
    add_value_matrix(solver, nconn, l_bi, l_bi, 2, 1, v_nconn)
    @. v_nconn = -grad(fw, 1, 2)
    add_value_matrix(solver, nconn, l_bi, l_bi, 1, 2, v_nconn)
    @. v_nconn = -grad(fo, 1, 2)
    add_value_matrix(solver, nconn, l_bi, l_bi, 2, 2, v_nconn)

    @. v_nconn = -grad(fw, 2, 1)
    add_value_matrix(solver, nconn, l_bi, r_bi, 1, 1, v_nconn)
    @. v_nconn = -grad(fo, 2, 1)
    add_value_matrix(solver, nconn, l_bi, r_bi, 2, 1, v_nconn)
    @. v_nconn = -grad(fw, 2, 2)
    add_value_matrix(solver, nconn, l_bi, r_bi, 1, 2, v_nconn)
    @. v_nconn = -grad(fo, 2, 2)
    add_value_matrix(solver, nconn, l_bi, r_bi, 2, 2, v_nconn)

    @. v_nconn = grad(fw, 1, 1)
    add_value_matrix(solver, nconn, r_bi, l_bi, 1, 1, v_nconn)
    @. v_nconn = grad(fo, 1, 1)
    add_value_matrix(solver, nconn, r_bi, l_bi, 2, 1, v_nconn)
    @. v_nconn = grad(fw, 1, 2)
    add_value_matrix(solver, nconn, r_bi, l_bi, 1, 2, v_nconn)
    @. v_nconn = grad(fo, 1, 2)
    add_value_matrix(solver, nconn, r_bi, l_bi, 2, 2, v_nconn)

    @. v_nconn = grad(fw, 2, 1)
    add_value_matrix(solver, nconn, r_bi, r_bi, 1, 1, v_nconn)
    @. v_nconn = grad(fo, 2, 1)
    add_value_matrix(solver, nconn, r_bi, r_bi, 2, 1, v_nconn)
    @. v_nconn = grad(fw, 2, 2)
    add_value_matrix(solver, nconn, r_bi, r_bi, 1, 2, v_nconn)
    @. v_nconn = grad(fo, 2, 2)
    add_value_matrix(solver, nconn, r_bi, r_bi, 2, 2, v_nconn)

    for well in values(wells)
        # Add well rate to residual
        ind = Int32.(well.ind)
        nperf = length(ind)
        ∂rw∂po, ∂rw∂sw, ∂ro∂po, ∂ro∂sw = -grad(well.qw, 1, 1), -grad(well.qw, 1, 2), -grad(well.qo, 1, 1), -grad(well.qo, 1, 2)
        
        # Add derivatives of reservoir equations to reservoir primary variables
        add_value_matrix(solver, nperf, ind, ind, 1, 1, ∂rw∂po) 
        add_value_matrix(solver, nperf, ind, ind, 1, 2, ∂rw∂sw) 
        add_value_matrix(solver, nperf, ind, ind, 2, 1, ∂ro∂po) 
        add_value_matrix(solver, nperf, ind, ind, 2, 2, ∂ro∂sw) 

        # Add derivatives of 
        # 1. reservoir equations to well variables
        # 2. well equations to reservoir primary variables
        # 3. well equations to well variables
        # 4. perform Schur complement in place : J -= Cw * Bw / Dw
        ∂rw∂w, ∂ro∂w = -grad(well.qw, 2, 1), -grad(well.qo, 2, 1) # Cw
        ∂w∂po, ∂w∂sw = well.∂r[:, 1], well.∂r[:, 2] # Bw
        ∂w∂w = well.∂w # Dw
        for i in ind
            BI = Int32(i) .* ones(Int32, nperf)
            add_value_matrix(solver, nperf, BI, ind, 1, 1, -∂rw∂w .* ∂w∂po ./ ∂w∂w) 
            add_value_matrix(solver, nperf, BI, ind, 1, 2, -∂rw∂w .* ∂w∂sw ./ ∂w∂w)  
            add_value_matrix(solver, nperf, BI, ind, 2, 1, -∂ro∂w .* ∂w∂po ./ ∂w∂w)  
            add_value_matrix(solver, nperf, BI, ind, 2, 2, -∂ro∂w .* ∂w∂sw ./ ∂w∂w) 
        end
    end
end  

function update_solution(nsolver::NRSolverDuneIstl, fluid::OWFluid, wells::Dict{String, AbstractFacility})
    o, w = fluid.phases.o, fluid.phases.w
    solver, assembler = nsolver.lsolver.solver, nsolver.assembler
    # get_value_x(int nn, int* BI, int I, T* value)
    get_value_x(solver, assembler.nc, assembler.diag_bi, 1, assembler.po)
    get_value_x(solver, assembler.nc, assembler.diag_bi, 2, assembler.sw)
    # Update well solutions
    for well in values(wells)
        ind = Int32.(well.ind)
        nperf = length(ind)
        ∂w∂po, ∂w∂sw = well.∂r[:, 1], well.∂r[:, 2] # Bw
        ∂w∂w = well.∂w # Dw
        rw = well.rw
        well.bhp -= (rw - sum(∂w∂po .* assembler.po[ind]) - sum(∂w∂sw .* assembler.sw[ind])) / ∂w∂w
    end

    # Apply appleyard chopping on saturation
    # println("Appleyard chopping")
    @. assembler.sw = min(assembler.sw, 0.35*value(w.s))

    @. assembler.po = value(o.p) - assembler.po
    @. assembler.sw = value(w.s) - assembler.sw
    # Bounds check
    @. assembler.sw = max(assembler.sw, 0.0)
    @. assembler.sw = min(assembler.sw, 1.0)
    update_primary_variable(fluid, assembler.po, assembler.sw)
    return fluid
end
