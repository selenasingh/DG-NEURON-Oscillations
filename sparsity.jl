include("utilities.jl");
default(show=false)
default(fontfamily="Computer Modern")

# HYPERPARAMS
n_runs = 3
patterns = 0:12
n_patterns = length(patterns)
duration = 2000
labels = ["theta" , "ftheta", "alpha", "beta", "gamma"]
freqlabels = [L"\theta" L"\theta_{fast}" L"\alpha" L"\beta"  L"\gamma"]

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

freqs = [L"\theta", L"\theta_{fast}", L"\alpha", L"\beta",  L"\gamma"]

spar_save = OrderedDict("theta"=>[], "ftheta"=>[], "alpha"=>[], "beta"=>[], "gamma"=>[])
active_gcs = OrderedDict("theta"=>[], "ftheta"=>[], "alpha"=>[], "beta"=>[], "gamma"=>[])

spar_means = OrderedDict("theta"=>[], "ftheta"=>[], "alpha"=>[], "beta"=>[], "gamma"=>[])
spar_ses = OrderedDict("theta"=>[], "ftheta"=>[], "alpha"=>[], "beta"=>[], "gamma"=>[])

activem = OrderedDict("theta"=>[], "ftheta"=>[], "alpha"=>[], "beta"=>[], "gamma"=>[])
actives = OrderedDict("theta"=>[], "ftheta"=>[], "alpha"=>[], "beta"=>[], "gamma"=>[])


for i ∈ 1:length(labels)
    for run ∈ 1:n_runs
        spikes = load_spike_files(patterns, labels[i]*"-$run", populations)
        for j ∈ 1:(n_patterns-1)
            gc_test = spikes[(spikes.Population .== "DG") .& (spikes.Pattern .== j), ["Neuron", "Time"]] 
            
            #Identify index for end of first theta cycle
            index = 1
            for i ∈ 1:length(gc_test[!, "Time"])
                timestamp = gc_test[i,"Time"]
                if timestamp < 285
                    index = index + 1
                end
            end

            #count up number of times a given GC is flagged for emitting a spike
            c = counter(gc_test[1:index,"Neuron"])
            vals = values(c) 
            sparsity = 1/mean(vals) 
            append!(spar_save[labels[i]], sparsity)

            popsparse = 500/length(vals)
            append!(active_gcs[labels[i]], popsparse)
        end
        if (run == n_runs) 
            sparm = round(mean(spar_save[labels[i]]), digits=2)
            append!(spar_means[labels[i]], sparm)
            sparse = std(spar_save[labels[i]])/sqrt(n_runs)
            append!(spar_ses[labels[i]], sparse)
            active_m = round(mean(active_gcs[labels[i]]), digits=2)
            append!(activem[labels[i]], active_m)
            active_s = std(active_gcs[labels[i]])/sqrt(n_runs)
            append!(actives[labels[i]], active_s)
        end
    end
end

unpack(a) = eltype(a[1])[el[1] for el in a]
sparsityGC_fig = plot(freqs, 
                unpack(collect(values(spar_means))), 
                xlabel = "Input Frequency Band",
                xtickfont=font(12),
                ylabel = "Sparsity of GC Firing",
                c = :black, 
                linewidth = 2,
                yerror = unpack(collect(values(spar_ses))), 
                dpi=300, size=(350,350),
                label=nothing,
                )

CSV.write("figures/pattern-separation/spar_means.csv", spar_means)
CSV.write("figures/pattern-separation/spar_ses.csv", spar_ses)

savefig(sparsityGC_fig, "figures/pattern-separation/sparsity-gc-curve"*fig_ext)
#sparsity_fig

sparsityNET_fig = plot(freqs, 
                unpack(collect(values(activem))), 
                xlabel = "Input Frequency Band",
                xtickfont=font(12),
                ylabel = "GC Population Sparsity",
                c = :black, 
                linewidth = 2,
                yerror = unpack(collect(values(actives))), 
                dpi=300, size=(350,350),
                label=nothing,
                )

CSV.write("figures/pattern-separation/net_spar_means.csv", activem)
CSV.write("figures/pattern-separation/net_spar_ses.csv", actives)
                
savefig(sparsityNET_fig, "figures/pattern-separation/sparsity-net-curve"*fig_ext)