module State

#! format: off
using ..Global: M
using ..AutoDiff: Tensor, param, zeros_tensor, data
using ..Fluid: AbstractPVTO, AbstractPVTW, AbstractSWOF, compute_kro,
        compute_krw, compute_bo, compute_bw, compute_μo, compute_μw
#! format: on
export OWState, AbstractState, set_init_state

abstract type AbstractState end

# Define struct State
struct OWState <: AbstractState
    numcell::Int
    numvar::Int
    ρo_std::Float64
    ρw_std::Float64
    p::Tensor # Pressure
    sw::Tensor  # Wate phase saturation

    so::Tensor # Oil phase saturation
    bo::Tensor # Oil phase formation volume factor
    bw::Tensor # Water phase formation volume factor
    μo::Tensor # Oil phase viscosity
    μw::Tensor # Water phase viscosity
    kro::Tensor # Oil phase rel-perm
    krw::Tensor # Water phase rel-perm
    λo::Tensor # Oil phase mobility
    λw::Tensor # Water phase mobility
    ρo::Tensor # Oil phase density
    ρw::Tensor # Water phase density
    ao::Tensor # Oil component accumulation
    aw::Tensor # Water component accumulation
    ro::Tensor # Oil component residual
    rw::Tensor # Water component residual

    γo::Tensor # Oil gravity term for each connection
    γw::Tensor # Oil gravity term for each connection
    ΔΨo::Tensor # Oil gravity potential difference for each connection
    ΔΨw::Tensor # Water gravity potential difference for each connection
    fo::Tensor # Oil component flux
    fw::Tensor # Water component flux

    pn::Vector{Float64}   # Pressure at previous time step
    swn::Vector{Float64} # Water saturation at last time step
    son::Vector{Float64} # Oil saturation at last time step
    bon::Vector{Float64} # bo at last time step
    bwn::Vector{Float64} # bw at last time step

    pvto::AbstractPVTO
    pvtw::AbstractPVTW
    swof::AbstractSWOF

    t::Vector{Float64}  # Save results
    p_rec::Vector{Vector{Float64}}
    sw_rec::Vector{Vector{Float64}}

    function OWState(
        nc::Int,
        nconn::Int,
        ρo_std::Float64,
        ρw_std::Float64,
        pvto::AbstractPVTO,
        pvtw::AbstractPVTW,
        swof::AbstractSWOF,
    )::OWState
        nv = 2
        p = param(zeros(nc), 1, nv)
        sw = param(zeros(nc), 2, nv)

        #! format: off
        vecs = [:so, :bo, :bw, :μo, :μw, :kro, :krw,
                    :λo, :λw, :ρo, :ρw, :ao, :aw, :ro, :rw]
        #! format: on
        tensors = [zeros_tensor(nc, nv) for v in vecs]

        vecs = [:γo, :γw, :ΔΨo, :ΔΨw, :fo, :fw]
        conn_tensors = [zeros_tensor(nconn, nv) for v in vecs]

        vecs = [:pn, :swn, :son, :bon, :bwn]
        params = [zeros(Float64, nc) for v in vecs]

        t = Vector{Float64}()
        p_rec = Vector{Vector{Float64}}()
        sw_rec = Vector{Vector{Float64}}()

        #! format: off
        return new(nc, nv, ρo_std, ρw_std, p, sw, tensors..., conn_tensors...,
                    params..., pvto, pvtw, swof, t, p_rec, sw_rec)
        #! format: on
    end
end

#
function get_var_order(state::OWState)::Dict{String, Int}
    return Dict{String, Int}("p" => 1, "sw" => 2)
end
#
function update_old_state(state::OWState)::OWState
    state.pn .= data(state.p)
    state.swn .= data(state.sw)
    state.son .= data(state.so)
    state.bon .= data(state.bo)
    state.bwn .= data(state.bw)
    return state
end
#
function save_state_result(state::OWState, t::Float64)::OWState
    push!(state.t, t)
    push!(state.p_rec, data(state.p))
    push!(state.sw_rec, data(state.sw))
    return state
end
#
function change_back_state(state::OWState)::OWState
    for i = 1:state.numcell
        state.p[i].val = state.pn[i]
        state.sw[i].val = state.swn[i]
    end
    return state
end

#
function set_init_state(state::OWState, p0::Float64, sw0::Float64)::OWState
    return set_init_state(state, p0*ones(state.numcell), sw0*ones(state.numcell))
end

function set_init_state(state::OWState, p0::Vector{Float64}, sw0::Vector{Float64})::OWState
    @assert all(p0 .> 0.0) && all(0.0 .<= sw0 .<= 1.0)
    state.p .= param(p0, 1, 2)
    state.sw .= param(sw0, 2, 2)
    compute_params(state)
    update_old_state(state)
end

function compute_so(so::Tensor, sw::Tensor)::Tensor
    so .= 1 .- sw
end

function compute_ρo(ρo::Tensor, bo::Tensor, ρo_std::Float64)::Tensor
    ρo .= ρo_std ./ bo
end

function compute_ρw(ρw::Tensor, bw::Tensor, ρw_std::Float64)::Tensor
    ρw .= ρw_std ./ bw
end

function compute_λo(λo::Tensor, kro::Tensor, μo::Tensor, bo::Tensor)::Tensor
    λo .= kro ./ (μo .* bo)
end

function compute_λw(λw::Tensor, krw::Tensor, μw::Tensor, bw::Tensor)::Tensor
    λw .= krw ./ (μw .* bw)
end



function compute_params(state::OWState)::Nothing
    compute_so(state.so, state.sw)
    compute_kro(state.swof, state.kro, state.so)
    compute_krw(state.swof, state.krw, state.sw)
    compute_bo(state.pvto, state.bo, state.p)
    compute_μo(state.pvto, state.μo, state.p)
    compute_bw(state.pvtw, state.bw, state.p)
    compute_μw(state.pvtw, state.μw, state.p)

    compute_ρo(state.ρo, state.bo, state.ρo_std)
    compute_ρw(state.ρw, state.bw, state.ρw_std)

    compute_λo(state.λo, state.kro, state.μo, state.bo)
    compute_λw(state.λw, state.krw, state.μw, state.bw)
    return nothing
end


end
