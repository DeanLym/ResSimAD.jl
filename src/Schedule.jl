module Schedule

using ..State:OWState
using ..AutoDiff:data

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
    t_end::Float64
    report_time::Vector{Float64}
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
        sch.t_end = Inf
        sch.report_time = Vector{Float64}()
        return sch
    end
end

function reset_time_step(sch::Scheduler)
    sch.dt = sch.dt0
    sch.dt_old = sch.dt
end

function update_dt(sch::Scheduler, state::OWState, converge::Bool)
    if !converge
        sch.dt *= 0.5
        sch.t_next = sch.t_current + sch.dt
        return nothing
    else
        ω, ηp, ηs = sch.ω, sch.ηp, sch.ηs
        sch.dt_old = sch.dt
        rp = (1+ω)*ηp ./ (abs.(data(state.p) .- state.pn) .+ ω*ηp)
        rs = (1+ω)*ηs ./ (abs.(data(state.so) .- state.son) .+ ω*ηs)
    end
    dt = min(sch.dt_old * min(minimum(rp), minimum(rs)), sch.dt_max)
    t_next = sch.t_next
    for t in sch.report_time
        if (t_next - t)*(t_next + dt - t) < 0
            dt = t - t_next
            break
        end
    end
    sch.dt = dt
    sch.t_current = sch.t_next
    sch.t_next = sch.t_current + dt
    return nothing
end

function set_report_time(scheduler::Scheduler, report_time::Vector{Float64})::Nothing
    scheduler.report_time = report_time
    return nothing
end

end
