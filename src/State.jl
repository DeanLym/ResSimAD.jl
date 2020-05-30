module State

using ..Global: M
using ..Grid: AbstractGrid
using ..AutoDiff: Tensor, param, zeros_tensor, data

export OWState, AbstractState, set_init_state

abstract type AbstractState end


# Define struct State
struct OWState <: AbstractState
    numcell::Int
    numvar::Int
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
    qo::Tensor # Oil component well rate
    qw::Tensor # Water component well rate
    ao::Tensor # Oil component accumulation
    aw::Tensor # Water component accumulation
    ro::Tensor # Oil component residual
    rw::Tensor # Water component residual

    fo::Tensor # Oil component flux
    fw::Tensor # Water component flux

    pn::Vector{Float64}   # Pressure at previous time step
    swn::Vector{Float64} # Water saturation at last time step
    son::Vector{Float64} # Oil saturation at last time step
    bon::Vector{Float64} # bo at last time step
    bwn::Vector{Float64} # bw at last time step

    function OWState(nc::Int, nconn::Int)::OWState
        nv = 2
        p = param(zeros(nc), 1, nv)
        sw = param(zeros(nc), 2, nv)

        #! format: off
        vecs = [:so, :bo, :bw, :μo, :μw, :kro, :krw,
                    :λo, :λw, :qo, :qw, :ao, :aw, :ro, :rw]
        #! format: on
        tensors = [zeros_tensor(nc, nv) for v in vecs]

        fo = zeros_tensor(nconn, nv) # Oil component flux
        fw = zeros_tensor(nconn, nv) # Water component flux

        vecs = [:pn, :swn, :son, :bon, :bwn]
        params = [zeros(Float64, nc) for v in vecs]

        #! format: off
        return new(nc, nv, p, sw, tensors..., fo, fw, params...)
        #! format: on
    end
end


#

function get_var_order(state::OWState)::Dict{String, Int}
    return Dict{String, Int}("p" => 1, "sw" => 2)
end
#
function update_old_state(state::OWState)::Nothing
    state.pn .= data(state.p)
    state.swn .= data(state.sw)
    state.son .= data(state.so)
    state.bon .= data(state.bo)
    state.bwn .= data(state.bw)
    return nothing
end

function change_back_state(state::OWState)::Nothing
    for i = 1:state.numcell
        state.p[i].val = state.pn[i]
        state.sw[i].val = state.swn[i]
    end
end

#
function set_init_state(state::OWState, p0::Float64, sw0::Float64)::Nothing
    return set_init_state(state, p0*ones(state.numcell), sw0*ones(state.numcell))
end

function set_init_state(state::OWState, p0::Vector{Float64}, sw0::Vector{Float64})::Nothing
    @assert all(p0 .> 0.0) && all(0.0 .<= sw0 .<= 1.0)
    state.p .= param(p0, 1, 2)
    state.sw .= param(sw0, 2, 2)
    compute_params(state)
    update_old_state(state)
    return nothing
end

function compute_so(so::Tensor, sw::Tensor)::Nothing
    so .= 1 .- sw
    return nothing
end

# # Define relative permeability curve
function compute_kro(kro::Tensor, so::Tensor)::Nothing
    kro .= so.^2
    return nothing
end
#
function compute_krw(krw::Tensor, sw::Tensor)::Nothing
    krw .= sw.^2
    return nothing
end
#
# # Define oil formation volume factor function
function compute_bo(bo::Tensor, po::Tensor)::Nothing
    co = 1.2e-5
    pref = 14.7
    bo .= @. exp(-co * (po - pref))
    return nothing
end
#
# # Define water formation volume factor function
function compute_bw(bw::Tensor, pw::Tensor)::Nothing
    cw = 5e-7
    pref = 14.7
    bw .= @. exp(-cw * (pw - pref))
    return nothing
end
#
# # Define viscosity function
function compute_μo(μo::Tensor, po::Tensor)::Nothing
    cμo = 2e-6
    pref = 14.7
    μo .= @. 5.0 * exp(cμo * (po - pref))
    return nothing
end
#
function compute_μw(μw::Tensor, po::Tensor)::Nothing
    μw .= 0 .* po .+ 1.0
    return nothing
end
#
function compute_λo(λo::Tensor, kro::Tensor, μo::Tensor, bo::Tensor)::Nothing
    λo .= @. kro / (μo * bo)
    return nothing
end
#
function compute_λw(λw::Tensor, krw::Tensor, μw::Tensor, bw::Tensor)::Nothing
    λw .= @. krw / (μw * bw)
    return nothing
end
#
#! format: off
function compute_ao(state::OWState, grid::AbstractGrid, dt::Float64)::Nothing
    state.ao .= 1.0 / M .* grid.v .* grid.ϕ .*
        (state.so ./ state.bo - state.son ./ state.bon) / dt
    return nothing
end
#! format: on

#! format: off
function compute_aw(state::OWState, grid::AbstractGrid, dt::Float64)::Nothing
    state.aw .= 1.0 / M .* grid.v .* grid.ϕ .*
        (state.sw ./ state.bw - state.swn ./ state.bwn) / dt
    return nothing
end
#! format: on

function compute_fo(state::OWState, grid::AbstractGrid)::Nothing
    trans = grid.connlist.trans
    for i = 1:grid.connlist.numconn
        l, r = grid.connlist.l[i], grid.connlist.r[i]
        if state.p[l] > state.p[r]
            state.fo[i] = state.λo[l] * trans[i] * (state.p[l] - state.p[r])
        else
            state.fo[i] = state.λo[r] * trans[i] * (state.p[l] - state.p[r])
        end
    end
    return nothing
end

function compute_fw(state::OWState, grid::AbstractGrid)::Nothing
    trans = grid.connlist.trans
    for i = 1:grid.connlist.numconn
        l, r = grid.connlist.l[i], grid.connlist.r[i]
        if state.p[l] > state.p[r]
            state.fw[i] = state.λw[l] * trans[i] * (state.p[l] - state.p[r])
        else
            state.fw[i] = state.λw[r] * trans[i] * (state.p[l] - state.p[r])
        end
    end
    return nothing
end

function compute_qo(state::OWState, prod_bhp::Vector{Float64})::Nothing
    wi = 0.5 # well index
    for i = 1:state.numcell
        if prod_bhp[i] != Inf
            state.qo[i] = wi * state.λo[i] * (state.p[i] - prod_bhp[i])
        end
    end
    return nothing
end

function compute_qw(state::OWState, prod_bhp::Vector{Float64}, inj_bhp::Vector{Float64})::Nothing
    wi = 0.5 # well index
    for i = 1:state.numcell
        if prod_bhp[i] != Inf
            state.qw[i] = wi * state.λw[i] * (state.p[i] - prod_bhp[i])
        end
        if inj_bhp[i] != Inf
            state.qw[i] = wi * state.λw[i] * (state.p[i] - inj_bhp[i])
        end
    end
    return nothing
end

function compute_params(state::OWState)::Nothing
    compute_so(state.so, state.sw)
    compute_kro(state.kro, state.so)
    compute_krw(state.krw, state.sw)
    compute_bo(state.bo, state.p)
    compute_bw(state.bw, state.p)
    compute_μo(state.μo, state.p)
    compute_μw(state.μw, state.p)

    compute_λo(state.λo, state.kro, state.μo, state.bo)
    compute_λw(state.λw, state.krw, state.μw, state.bw)
    return nothing
end


function compute_ro(
    state::OWState,
    grid::AbstractGrid,
    prod_bhp::Vector{Float64},
    dt::Float64,
)::Nothing
    compute_ao(state, grid, dt)
    compute_fo(state, grid)
    compute_qo(state, prod_bhp)
    state.ro .= - state.ao .- state.qo
    l, r = grid.connlist.l, grid.connlist.r
    for i = 1:grid.connlist.numconn
        state.ro[l[i]] -= state.fo[i]
        state.ro[r[i]] += state.fo[i]
    end
    return nothing
end

function compute_rw(
    state::OWState,
    grid::AbstractGrid,
    prod_bhp::Vector{Float64},
    inj_bhp::Vector{Float64},
    dt::Float64,
)::Nothing
    compute_aw(state, grid, dt)
    compute_fw(state, grid)
    compute_qw(state, prod_bhp, inj_bhp)
    state.rw .= - state.aw .- state.qw
    l, r = grid.connlist.l, grid.connlist.r
    for i = 1:grid.connlist.numconn
        state.rw[l[i]] -= state.fw[i]
        state.rw[r[i]] += state.fw[i]
    end
    return nothing
end

end
