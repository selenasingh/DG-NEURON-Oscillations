## LFP CALCULATION AND PLOTTING POWER SPECTRA
## - built from referencing 'spike_train_correlation' in utilities.jl
## ---------------------------------------------

include("utilities.jl");
default(show=false)
default(fontfamily="Computer Modern")

# HYPERPARAMS
n_runs = 2
patterns = 0:12
duration = 750
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


# FUNCTIONS

function makechunks(convolution::Vector, n_neurons::Integer)
    chunks = length(convolution) ÷ n_neurons
    return [convolution[1+chunks*k:(k == n_neurons-1 ? end : chunks*k+chunks)] for k = 0:n_neurons-1]
end  

function compute_LFP(
    spike_train::DataFrame, 
    n_neurons::Int64;
    delta::Int64=10, 
    n_bins::Int64=1,        
    duration::Float64=750.0)
    
    # Create kernel
    tri_rng = collect(1:round(Int64, delta/n_bins))
    triangular_kernel = [tri_rng; reverse(tri_rng)[2:end]]'

    # Create bins for histograms 
    bins = collect(1:n_bins:(duration + n_bins))

    convolved_series = zeros(round(Int64, n_neurons * duration / n_bins)) 
    for i ∈ 1:n_neurons
        s = spike_train[spike_train.Neuron .== (i-1), "Time"]
        timeseries, _ = np.histogram(s, bins)
        idx_low = round(Int64, (i-1)*duration/n_bins + 1)
        idx_high = round(Int64, i*duration/n_bins)
        convolved_series[idx_low:idx_high] = np.convolve(timeseries, triangular_kernel, "same")     #TODO: consider scipy.signal.fftconvolve as alt.
    end
    LFP = makechunks(convolved_series, n_neurons)
    return sum(LFP)
end

global ppsave = []
global gcsave = []
lfp_save_pp = OrderedDict("theta"=>[], "ftheta"=>[], "alpha"=>[], "beta"=>[], "gamma"=>[])
lfp_save_gc = OrderedDict("theta"=>[], "ftheta"=>[], "alpha"=>[], "beta"=>[], "gamma"=>[])
for run_ ∈ 1:n_runs
    for i ∈ 1:length(labels)
        spikes = load_spike_files(patterns, labels[i]*"-$run_", populations)

        for p ∈ unique(spikes.Pattern)
            # TODO: change so this works ∀ neurons
            pp = spikes[(spikes.Population .== "PP") .& (spikes.Pattern .== p), :]
            gc = spikes[(spikes.Population .== "DG") .& (spikes.Pattern .== p), :]
            
            #compute LFP
            pp_lfp = [compute_LFP(pp, 100)]
            gc_lfp = [compute_LFP(gc, 500)]

            #save LFP for each pattern:
            global ppsave = append!(ppsave, pp_lfp)
            global gcsave = append!(gcsave, gc_lfp)
        end
        # average LFP across all patterns
        sum_lfp_pat_pp = [sum(ppsave)/length(ppsave[1])]
        sum_lfp_pat_gc = [sum(gcsave)/length(gcsave[1])]

        # save summated lfp/freq band
        append!(lfp_save_pp[labels[i]], sum_lfp_pat_pp)
        append!(lfp_save_gc[labels[i]], sum_lfp_pat_gc)
    end
end

l = @layout [grid(2,1) a{0.5w}]
theta_lfp_pp = plot(lfp_save_pp["theta"][1], 
            color = :black,
            xlabel = "Time (ms)",
            ylabel = "Population Activity",
            label = L"\theta_{PP}"
)          
#savefig(fig, "figures/lfp/test_pp.png")

theta_lfp_gc = plot(lfp_save_gc["theta"][1], 
            color = :black,
            xlabel = "Time (ms)",
            ylabel = "Population Activity",
            label = L"GC"
)
#savefig(fig, "figures/lfp/test_gc.png")

# POWER SPECTRA 
restruct_gclfps = zeros((duration, length(labels)))
restruct_pplfps = zeros((duration, length(labels)))

for i ∈ 1:length(labels)
    restruct_gclfps[:,i] = lfp_save_gc[labels[i]][n_runs]
    restruct_pplfps[:,i] = lfp_save_pp[labels[i]][n_runs]
end

S = spectra(restruct_gclfps, duration, duration; tapering=hamming) #TODO: learn what this 'tapering' is.
powerfig = plot(S, fmax = 40, label = freqlabels, legend = :topright, ylabel = "Power "*L"(\mu V^2)")
#savefig(powerfig, "figures/lfp/powerspectra_gcs.png")

sumfig = plot(theta_lfp_pp, theta_lfp_gc, powerfig, layout=l)
savefig(sumfig, "figures/lfp/powerspectra_lfps_gcs.png")


