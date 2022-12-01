## WORKING TOWARDS LFP CALCULATION 
## - built from referencing 'spike_train_correlation' in utilities.jl
## - Need to figure out where to truncate convolution 
## ---------------------------------------------

include("utilities.jl");
default(show=false)

# HYPERPARAMS
n_runs = 1
patterns = 0:1
labels = ["theta", "alpha", "gamma"]
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

function compute_LFP(
    spike_train::DataFrame, 
    n_neurons::Int64;
    delta::Int64=100, 
    n_bins::Int64=1,        #TODO: understand why there's only 1 bin..? What's the relationship with delta?
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
    return convolved_series 
end

# Load theta file only
spikes = load_spike_files(patterns, labels[1]*"-1", populations)

# separate above df into neuronal populations
pp_a = spikes[(spikes.Population .== "PP") .& (spikes.Pattern .== 1), ["Neuron", "Time"]]
gc_test = spikes[(spikes.Population .== "DG") .& (spikes.Pattern .== 1), ["Neuron", "Time"]]

# testing compute_LFP f'n above
pp_lfp = compute_LFP(pp_a, 100)
fig = plot(pp_lfp[1:750])       # this looks v promising
savefig(fig, "figures/lfp/test_pp.png")

gc_lfp = compute_LFP(gc_test, 500)
convolution_size=length(gc_lfp)
midpoint = div(convolution_size, 2)     #returns int

#print(convolution_size)
#print(midpoint)

#TODO: understand where would be best to truncate, how to interpret axes
fig = plot(gc_lfp[(midpoint-7500):(midpoint+7500)])   
savefig(fig, "figures/lfp/test_gc.png")

