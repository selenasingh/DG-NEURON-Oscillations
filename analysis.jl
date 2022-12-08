
####################################################################
# SCRIPT TO RUN ANALYSIS 
####################################################################
include("utilities.jl");
default(show=false)

# HYPERPARAMS
n_runs = 2
patterns = 0:12
labels = ["theta" , "ftheta", "alpha", "beta", "gamma"]
freqs = [L"\theta", L"\theta_{fast}", L"\alpha", L"\beta",  L"\gamma"]
fig_ext = ".png"

default(fontfamily="Computer Modern")

# CREATE NECESSARY DIRECTORIES 
create_directories_if_not_exist()

# IDENTIFY NEURON POPULATION RANGES
populations = Dict(
    "DG" => [0, 500],
    "BC" => [500,506],
    "MC" => [506, 521],
    "HIPP" => [521, 527]
)


for run_ ∈ 1:n_runs
    for i ∈ 1:length(labels)
        spikes = load_spike_files(patterns, labels[i]*"-$run_", populations)

        # CREATE RASTER PLOTS
        for p ∈ unique(spikes.Pattern)
            stimin = spikes[(spikes.Population .== "PP") .& (spikes.Pattern .== p), :]
            plots = []
            append!(plots, [raster_plot(stimin; ylab="PP")])
            for pop ∈ keys(populations)
                lb, ub = populations[pop]
                popspikes = spikes[(spikes.Population .== pop) .& (spikes.Pattern .== p),:]
                append!(plots, [raster_plot(popspikes; xlab="", ylab=pop)])
            end
            fig = plot(reverse(plots)..., layout=grid(5, 1, heights=[0.15, 0.15, 0.15, 0.4, 0.15]), size=(400, 500))
            savefig(fig, "figures/raster-plots/raster-"*string(p)*"-"*labels[i]*"-$run_"*fig_ext)
        end
    end 
end 


# PATTERN SEPARATION CURVES
colors=[:blue, :red, :green, :grey, :black]
global psfig = plot([0;1], [0;1], ls=:dash, c=:black, 
                        xlabel="Input Correlation "*L"(r_{in})", 
                        ylabel="Output Correlation "*L"(r_{out})", 
                        dpi=300, size=(350,350),  # increase resolution of image
                        label=nothing, legend=:topleft)

#psc = Dict("theta"=>[], "alpha"=>[], "gamma"=>[])
psc = Dict("theta"=>[], "ftheta"=>[], "alpha"=>[], "beta"=>[], "gamma"=>[])
for i ∈ 1:length(labels)
    for run ∈ 1:n_runs
        spikes = load_spike_files(patterns, labels[i]*"-$run", populations)
        
        out = pattern_separation_curve(spikes, 100, 500)
        x, y = out[:,"Input Correlation"], out[:, "Output Correlation"]
        
        # Remove NaNs before fitting
        idx_ = (.!isnan.(x) .& .!isnan.(y))
        x = x[idx_]
        y = y[idx_]

        f = fit_power_law(x, y)
        append!(psc[labels[i]], f(0.6))

        if (run == n_runs) 
            psm = round(mean(psc[labels[i]]), digits=2)
            psse = std(psc[labels[i]])/sqrt(n_runs)
            pslci = round(psm - 1.96*psse, digits=2)
            psuci = round(psm + 1.96*psse, digits=2)
            if n_runs > 1
                psc_label = freqs[i]*" (PS="*string(psm)*" ["*string(pslci)*", "*string(psuci)*"])"
            else 
                psc_label = freqs[i]*" (PS="*string(psm)*")"
            end
        else
            psc_label = nothing
        end
        global psfig = scatter!(x, y, c=colors[i], alpha=1/(2*n_runs), label=nothing)
        global psfig = plot!(0:0.01:1, x -> f(x), c=colors[i], label=psc_label)
    end 
end 
psfig
savefig(psfig, "figures/pattern-separation/pattern-separation-curve"*fig_ext)

# AREA UNDER PS CURVES 

#=
auc_save = OrderedDict("theta"=>[], "alpha"=>[], "gamma"=>[])
auc_means = OrderedDict("theta"=>[], "alpha"=>[], "gamma"=>[])
auc_ses = OrderedDict("theta"=>[], "alpha"=>[], "gamma"=>[])
=#

auc_save = OrderedDict("theta"=>[], "ftheta"=>[], "alpha"=>[], "beta"=>[], "gamma"=>[])
auc_means = OrderedDict("theta"=>[], "ftheta"=>[], "alpha"=>[], "beta"=>[], "gamma"=>[])
auc_ses = OrderedDict("theta"=>[], "ftheta"=>[], "alpha"=>[], "beta"=>[], "gamma"=>[])

for i ∈ 1:length(labels)
    for run ∈ 1:n_runs
        spikes = load_spike_files(patterns, labels[i]*"-$run", populations)
        
        out = pattern_separation_curve(spikes, 100, 500)
        x, y = out[:,"Input Correlation"], out[:, "Output Correlation"]
        
        # Remove NaNs before fitting
        idx_ = (.!isnan.(x) .& .!isnan.(y))
        x = x[idx_]
        y = y[idx_]

        auc = compute_auc(x, y)
        append!(auc_save[labels[i]], auc)
        if (run == n_runs) 
            aucm = round(mean(auc_save[labels[i]]), digits=2)
            append!(auc_means[labels[i]], aucm)
            aucse = std(auc_save[labels[i]])/sqrt(n_runs)
            append!(auc_ses[labels[i]], aucse)
        end
    end 
end 

unpack(a) = eltype(a[1])[el[1] for el in a]
auc_fig = plot(freqs, 
                unpack(collect(values(auc_means))), 
                xlabel = "Frequency Band",
                xtickfont=font(12),
                ylabel = L"AUC_{PS}",
                c = :black, 
                linewidth = 2,
                yerror = unpack(collect(values(auc_ses))), 
                dpi=300, size=(300,300),
                label=nothing,
                )
savefig(auc_fig, "figures/pattern-separation/auc-curve"*fig_ext)
