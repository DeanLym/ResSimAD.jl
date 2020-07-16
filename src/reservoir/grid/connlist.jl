mutable struct ConnList
    nconn::Int
    l::Vector{Int}
    r::Vector{Int}
    trans::Vector{Float64}
    Δd::Vector{Float64} # Depth Diffrence
    function ConnList()
        l = Int[]
        r = Int[]
        trans = Float64[]
        Δd = Float64[]
        return new(0, l, r, trans, Δd)
    end
end

function add_conn(conn::ConnList, l::Int, r::Int, trans::Float64, Δd::Float64)::ConnList
    push!(conn.l, l)
    push!(conn.r, r)
    push!(conn.trans, trans)
    push!(conn.Δd, Δd)
    return conn
end

function sort_conn(conn::ConnList)
    # Force l[i] < r[i]
    l, r = conn.l, conn.r
    for i in eachindex(l)
        if r[i] < l[i]
            l[i], r[i] = r[i], l[i]
        end
    end
end
