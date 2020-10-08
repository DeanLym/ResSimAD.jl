## Define some function for convenience

# Define function po(sim), pw(sim), so(sim), sw(sim) ....
for v in (:p, :s, :b, :μ, :kr, :λ, :ρ, :γ, :ΔΨ, :f, :pn, :sn, :bn, :ρs, :pvt)
    for p in ("o", "w", "g")
        @eval function $(Symbol(v,p))(x::Sim)
            getfield(x, :reservoir).fluid.phases[Symbol($p)].$v
        end
    end
end

# Define function ao(sim), aw(sim), ro(sim), rw(sim) ....
for v in (:a, :r)
    for p in ("o", "w", "g")
        @eval function $(Symbol(v,p))(x::Sim)
            getfield(x, :reservoir).fluid.components[Symbol($p)].$v
        end
    end
end

# Define function po_rec(sim), pw_rec(sim) ...
for v in ("p", "s")
    for p in ("o", "w", "g")
        @eval function $(Symbol(v, p, "_rec"))(x::Sim)
            getfield(getfield(x, :reservoir).fluid.phases[Symbol($p)], Symbol($v, "_rec"))
        end
    end
end

# Define function kx(sim), ϕ(sim) ...
for v in (:kx, :ky, :kz, :ϕ)
    @eval function $v(x::Sim)
        getfield(x, :reservoir).rock.$v
    end
end

# Define function nc(sim), nx(sim) ...
for v in (:nc, :nx, :ny, :nz, :dx, :dy, :dz, :v, :d, :connlist)
    @eval function $v(x::Sim)
        getfield(x, :reservoir).grid.$v
    end
end

# Define function residual, jacobian
for v in (:residual, :jac)
    @eval function $v(x::Sim)
        getfield(x, :nsolver).$v
    end
end

# Define dt(sim), t_current(sim), t_next(sim)
for v in (:dt, :t_current, :t_next)
    @eval function $v(x::Sim)
        getfield(x, :scheduler).$v
    end
end

# Define function reservoir(sim), facility(sim) ...
for v in ("reservoir", "facility", "scheduler", "nsolver")
    @eval $(Symbol(v))(x::Sim) = getfield(x, Symbol($v))
end
# Overwrite getproperty function for sim
# sim.$symbol = symbol(sim)
# Now we can use sim.po to get sim.reservoir.fluid.phases.o.p
# Similarly sim.ro => sim.reservoir.fluid.components.o.r
# sim.po_rec => sim.reservoir.fluid.phases.o.p_rec
# At the same time, sim.reservoir sim.facility sim.nsolver ... still works

import Base
Base.getproperty(x::Sim, s::Symbol) = eval(s)(x)
