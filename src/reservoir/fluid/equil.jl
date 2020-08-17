function cal_pressure(p_ref, δd, ρs, pvt, max_err)
    ρ1 = ρs/pvt.b(p_ref)
    p_old = p_ref + ρ1*δd*β
    ρ2 = ρs/pvt.b(p_old)
    ρ_avg = (ρ1 + ρ2) / 2.
    p_new = p_ref + ρ_avg*δd*β
    while abs(p_old - p_new) > max_err
        p_old = p_new
        ρ_avg = (ρs/pvt.b(p_ref)  + ρs/pvt.b(p_new)) / 2.
        # ρ_avg = ρs / pvt.b((p_ref + p_new) / 2.)
        p_new = p_ref + ρ_avg*δd*β
    end
    return p_new
end

function cal_pressure_table(num_node, d_ref, p_ref, δd, ρs, pvt, max_err)
    if num_node > 1
        d_table = Vector{Float64}(undef, num_node)
        for i = 1:num_node
            d_table[i] = d_ref + δd*(i-1)
        end
        p_new = cal_pressure(p_ref, δd, ρs, pvt, max_err)
        p_table = Vector{Float64}(undef, num_node)
        p_table[1] = p_ref
        p_table[2] = p_new
        for i = 3:num_node
            p_table[i] = cal_pressure(p_table[i-1], δd, ρs, pvt, max_err)
        end
    else
        d_table = [d_ref,]
        p_table = [p_ref,]
    end
    return d_table, p_table
end

function equil(
    fluid::AbstractFluid,
    d_ref::Float64,
    p_ref::Float64,
    d::Vector{Float64},
    Δd::Float64,
)
    d_min, d_max = minimum(d), maximum(d)

    max_err = 1.0e-6

    ρs = fluid.phases.o.ρs

    pvt = fluid.phases.o.pvt

    num_node_down = Int(ceil((d_max - d_ref) / Δd) + 1)
    num_node_up = Int(ceil((d_ref - d_min) / Δd) + 1)

    d_down, p_down = cal_pressure_table(num_node_down, d_ref, p_ref, Δd, ρs, pvt, max_err)
    d_up, p_up = cal_pressure_table(num_node_up, d_ref, p_ref, -Δd, ρs, pvt, max_err)

    d_table = vcat(reverse(d_up), d_down[2:end])
    p_table = vcat(reverse(p_up), p_down[2:end])

    p_interp = interpolate((d_table,), p_table, Gridded(Linear()))
    p_interp = extrapolate(p_interp, Flat())
    return p_interp.(d)
end
