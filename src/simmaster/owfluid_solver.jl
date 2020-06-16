
function compute_residual(
    fluid::OWFluid,
    grid::AbstractGrid,
    rock::AbstractRock,
    wells::Dict{String, AbstractFacility},
    dt::Float64
)
    oc, wc = fluid.components.o, fluid.components.w
    o, w = fluid.phases.o, fluid.phases.w
    ro, rw = oc.r, wc.r
    # Accumulation term
    ao, aw = oc.a, wc.a
    v, ϕ = grid.v, rock.ϕ
    ao .= compute_a(o, v, ϕ, dt)
    aw .= compute_a(w, v, ϕ, dt)

    ro .= -ao
    rw .= -aw

    # Well Term
    for well in values(wells)
        compute_qo(well, fluid)
        compute_qw(well, fluid)
        # Add well rate to residual
        ind = well.ind
        ro[ind] .-= well.qo
        rw[ind] .-= well.qw
    end

    # Flux term
    connlist = grid.connlist
    l, r = connlist.l, connlist.r
    nc, nv = grid.nc, fluid.nv
    fo, fw = o.f, w.f
    for i = 1:connlist.nconn
        ro[l[i]] -= fo[i]
        ro[r[i]] += fo[i]
        rw[l[i]] -= fw[i]
        rw[r[i]] += fw[i]
    end
    return fluid
end


function compute_residual_error(fluid::OWFluid, grid::AbstractGrid, rock::AbstractRock, dt::Float64)
    a = dt .* M ./ (grid.v .* rock.ϕ)
    op, wp = fluid.phases.o, fluid.phases.w
    oc, wc = fluid.components.o, fluid.components.w
    rw_err = a .* data(wp.b) .* data(wc.r)
    ro_err = a .* data(op.b) .* data(oc.r)
    return max(maximum(abs.(rw_err)), maximum((abs.(ro_err))))
end

function assemble_residual(fluid::OWFluid)::Vector{Float64}
    o, w = fluid.components.o, fluid.components.w
    ro, rw = o.r, w.r
    nc = length(ro)
    residual = zeros(Float64, 2*nc)

    for i = 1:nc
        residual[2*i-1] = rw[i].val
        residual[2*i] = ro[i].val
    end
    return residual
end

function assemble_jacobian(fluid::OWFluid)::SparseMatrixCSC{Float64,Int}
    I, J, V = Int[], Int[], Float64[]
    nv = fluid.nv
    o, w = fluid.components.o, fluid.components.w
    rw, ro = w.r, o.r
    nc = length(ro)
    for i =1:nc
        (ind, g) = grad(rw[i])
        n = length(ind)
        v1 = ones(Int, 2*n)
        v2 = ones(n)
        v3 = nv*ind
        append!(I, (2*i-1)*v1)
        for j = 1:nv
            append!(J, v3 + (j-nv)*v2)
        end
        append!(V, collect(transpose(g).data))

        (ind, g) = grad(ro[i])
        v3 = nv*ind
        append!(I, 2*i*v1)
        for j = 1:nv
            append!(J,v3 + (j-nv)*v2)
        end
        append!(V, collect(transpose(g).data))
    end
    return sparse(I, J, V);
end

function update_solution(fluid::OWFluid, δx::Vector{Float64})
    o, w = fluid.phases.o, fluid.phases.w
    po, sw = data(o.p), data(w.s)
    nc, nv = length(po), fluid.nv
    for i=1:nc
        po[i] -= δx[nv*(i-1) + 1]
        sw[i] -= δx[nv*(i-1) + 2]
    end
    # Bounds check
    for i =1:nc
        if sw[i] < 0.0
            sw[i] = 0.0
        elseif sw[i] > 1.0
            sw[i] = 1.0
        end
    end
    update_primary_variable(fluid, po, sw)
    return fluid
end
