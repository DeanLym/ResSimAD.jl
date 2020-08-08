using ResSimAD: get_grid_index
using ResSimAD: isproducer
## Define function to plot wells

function plot_wells_2d(plt, sim; markersize=5, color=:white)
    for w in values(sim.facility)
        i, j, _ = get_grid_index(sim.reservoir.grid, w.ind[1])
        if isproducer(w)
            marker = :circle
        else
            marker = :dtriangle
        end
        scatter!(plt, [j,], [i,], m=(marker, markersize, color), legend=false)
    end
end
