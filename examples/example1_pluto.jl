### A Pluto.jl notebook ###
# v0.11.6

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ dcd7cd30-e0ee-11ea-1c2b-b9915a81db2b
begin
	using ResSimAD
	using Plots
	using Plots.PlotMeasures
	include("utils.jl");
	cmap = cgrad(:jet);
	using PlutoUI
	md"##### Import modules"
end

# ╔═╡ 7b4cd170-e0f6-11ea-324d-79ae5b0779c3
md"###### This is a Pluto.jl notebook demonstrating interactive visualization of ResSimAD results."

# ╔═╡ dcfc2140-e0f2-11ea-2be0-ffe5feb97a19
md"##### Run simulation"

# ╔═╡ c51fb810-e0ee-11ea-3510-5b56410a99eb
sim, options = get_model("example1");

# ╔═╡ 5d1d6ff0-e0ee-11ea-0877-2fa6aa73502b
runsim(sim)

# ╔═╡ f5bcdeb0-e0f0-11ea-0dc6-51c2b34c794b
md"###### Plot simulation results"

# ╔═╡ 0f01ae50-e0f1-11ea-2490-6958ff5cf91d
begin
	wnames_ = collect(keys(sim.facility))
	h1 = @bind well_names PlutoUI.MultiSelect(wnames_; default=[wnames_[1]])
	h2 = @bind data_types PlutoUI.MultiSelect(["ORAT", "WRAT"], default=["WRAT"])
	md"""
	**Well names:** $(h1) --|-- **Data types**: $(h2)"""
end

# ╔═╡ 78d2ca60-e0ee-11ea-1902-49a8d8ee2166
begin
	t = get_well_rates(sim, well_names[1], "TIME");
	plt_handle = plot(size=(460, 280), legend=:outerright, xlabel="Day")
	for w in well_names
		for data in data_types
			y = get_well_rates(sim, w, data)
			plot!(plt_handle, t, abs.(y), label="$w $data", marker=true)
		end
	end
	plt_handle
end

# ╔═╡ c5134a10-e0f0-11ea-3d82-919c7b3b3957
md"###### Plot state maps"

# ╔═╡ ca71cd10-e0f0-11ea-2fe0-f17adbe47934
begin
	hslider = @bind tstep Slider(1:length(t))
	md"""
	Time step $(hslider)"""
end

# ╔═╡ 94c40a6e-e0f0-11ea-3951-b1d30a2adff7
begin
	# Plot state maps
	po = get_state_map(sim, "po", t[tstep]);
	sw = get_state_map(sim, "sw", t[tstep]);
	p3 = heatmap(reshape(po, sim.nx, sim.ny), color=cmap, title="Po at day $(round(t[tstep], digits=3))", clim=(5500., 6500.0));
	plot_wells_2d(p3, sim);
	p4 = heatmap(reshape(sw, sim.nx, sim.ny), color=cmap, title="Sw at day $(round(t[tstep], digits=3))", clim=(0., 1.0));
	plot_wells_2d(p4, sim);
	plot(p3, p4, layout=(1,2), size=(500,350))
end

# ╔═╡ Cell order:
# ╟─7b4cd170-e0f6-11ea-324d-79ae5b0779c3
# ╠═dcd7cd30-e0ee-11ea-1c2b-b9915a81db2b
# ╟─dcfc2140-e0f2-11ea-2be0-ffe5feb97a19
# ╠═c51fb810-e0ee-11ea-3510-5b56410a99eb
# ╠═5d1d6ff0-e0ee-11ea-0877-2fa6aa73502b
# ╟─f5bcdeb0-e0f0-11ea-0dc6-51c2b34c794b
# ╟─0f01ae50-e0f1-11ea-2490-6958ff5cf91d
# ╟─78d2ca60-e0ee-11ea-1902-49a8d8ee2166
# ╟─c5134a10-e0f0-11ea-3d82-919c7b3b3957
# ╟─ca71cd10-e0f0-11ea-2fe0-f17adbe47934
# ╟─94c40a6e-e0f0-11ea-3951-b1d30a2adff7
