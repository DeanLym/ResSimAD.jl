module Schedule

using ..AutoDiff:value
using ..Fluid:OWFluid, SPFluid


mutable struct Scheduler
    t0::Float64
    dt_old::Float64
    dt::Float64
    dt_max::Float64
    dt0::Float64
    ω::Float64
    ηp::Float64
    ηs::Float64
    t_current::Float64
    t_next::Float64
    time_step::Vector{Float64}
    function Scheduler()
        sch = new()
        sch.t0 = 0.0
        sch.dt0 = 0.01
        sch.dt = sch.dt0
        sch.dt_old = sch.dt
        sch.dt_max = Inf
        sch.ω = 0.5
        sch.ηp = 1000.
        sch.ηs = 0.4
        sch.t_current = 0.0
        sch.t_next = sch.t_current + sch.dt
        sch.time_step = [100.]
        return sch
    end
end

function reset_dt(sch::Scheduler)
    sch.dt = sch.dt0
    sch.dt_old = sch.dt
end

function update_dt(sch::Scheduler, fluid::OWFluid, converge::Bool)
    if !converge
        sch.dt *= 0.5
        sch.t_next = sch.t_current + sch.dt
        return nothing
    else
        ω, ηp, ηs = sch.ω, sch.ηp, sch.ηs
        sch.dt_old = sch.dt
        o = fluid.phases.o
        rp = (1+ω)*ηp ./ (abs.(value(o.p) .- o.pn) .+ ω*ηp)
        rs = (1+ω)*ηs ./ (abs.(value(o.s) .- o.sn) .+ ω*ηs)
    end
    dt = min(sch.dt_old * min(minimum(rp), minimum(rs)), sch.dt_max)
    t_next = sch.t_next
    for t in sch.time_step
        if (t_next - t)*(t_next + dt - t) < 0
            dt = t - t_next
            break
        end
    end
    sch.dt = dt
    sch.t_current = sch.t_next
    sch.t_next = round(sch.t_current + dt, digits=10)
    return nothing
end

function set_time_step(scheduler::Scheduler, time_step::Vector{Float64})::Vector{Float64}
    scheduler.time_step = time_step
end

function insert_time_step(scheduler::Scheduler, t::Float64)::Vector{Float64}
    push!(scheduler.time_step, t)
    sort!(scheduler.time_step)
end


function update_dt(sch::Scheduler, fluid::SPFluid, converge::Bool)
    if !converge
        sch.dt *= 0.5
        sch.t_next = sch.t_current + sch.dt
        return nothing
    else
        ω, ηp = sch.ω, sch.ηp
        sch.dt_old = sch.dt
        phase = fluid.phases[1]
        rp = minimum((1+ω)*ηp ./ (abs.(value(phase.p) .- phase.pn) .+ ω*ηp))
    end
    dt = min(sch.dt_old * rp, sch.dt_max)
    t_next = sch.t_next
    for t in sch.time_step
        if (t_next - t)*(t_next + dt - t) < 0
            dt = t - t_next
            break
        end
    end
    sch.dt = dt
    sch.t_current = sch.t_next
    sch.t_next = round(sch.t_current + dt, digits=10)
    return nothing
end

end
