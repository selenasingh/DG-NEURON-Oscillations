
####################################################################
# SCRIPT TO RUN ANALYSIS 
####################################################################
include("utilities.jl");


# HYPERPARAMS
patterns = 0:12
model_label = "test"
fig_ext = ".png"

# CREATE NECESSARY DIRECTORIES 
create_directories_if_not_exist()

# IDENTIFY NEURON POPULATION RANGES
populations = Dict(
    "DG" => [0, 500],
    "BC" => [500,506],
    "MC" => [506, 521],
    "HIPP" => [521, 527]
)

# LOAD DATA 
spikes = load_spike_files(patterns, model_label, populations)
spikes_spr = load_spike_files(patterns, "sprouted", populations)

# CREATE RASTER PLOTS
for p ∈ unique(spikes.Pattern)
    stimin = spikes[(spikes.Population .== "PP") .& (spikes.Pattern .== p), :]
    plots = []
    append!(plots, [raster_plot(stimin; ylab="PP")])
    for pop ∈ keys(populations)
        lb, ub = populations[pop]
        popspikes = spikes[(spikes.Population .== pop) .& (spikes.Pattern .== p),:]
        #if size(popspikes,1) > 0
        append!(plots, [raster_plot(popspikes; xlab="", ylab=pop)])
        #end
    end
    fig = plot(reverse(plots)..., layout=grid(5, 1, heights=[0.15, 0.15, 0.15, 0.4, 0.15]), size=(400, 500))
    savefig(fig, "figures/raster-"*string(p)*"-"*model_label*fig_ext)
end


out = pattern_separation_curve(spikes, 100, 500)
out_spr = pattern_separation_curve(spikes_spr, 100, 500)

x, y = out[:,"Input Correlation"], out[:, "Output Correlation"]
xs,ys = out_spr[:,"Input Correlation"], out_spr[:, "Output Correlation"]



f = fit_power_law(x, y)
f_spr = fit_power_law(xs, ys)
psc = round(f(0.6), digits=2)
psc2 = round(f_spr(0.6), digits=2)
round(0.555, digits=2)
using LaTeXStrings
psfig = plot([0;1], [0;1], ls=:dash, c=:black, xlabel=L"R_{in}", ylabel=L"R_{out}", label=nothing, legend=:bottomright)
psfig = scatter!(x, y, c=:blue, label="Base (PS="*string(psc)*")")
psfig = scatter!(xs, ys, c=:red, label="MFS 10% (PS="*string(psc2)*")")
psfig = plot!(0:0.01:1, x -> f(x), c=:blue, label=nothing)
psfig = plot!(0:0.01:1, x -> f_spr(x), c=:red, label=nothing)
