## LFP CALCULATION AND PLOTTING POWER SPECTRA
## - built from referencing 'spike_train_correlation' in utilities.jl
## ---------------------------------------------

include("utilities.jl");
default(show=false)
default(fontfamily="Computer Modern")

# HYPERPARAMS
n_runs = 3
patterns = 0:12
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


# FUNCTIONS

function makechunks(convolution::Vector, n_neurons::Integer)
    chunks = length(convolution) ÷ n_neurons
    return [convolution[1+chunks*k:(k == n_neurons-1 ? end : chunks*k+chunks)] for k = 0:n_neurons-1]
end  

function compute_LFP(
    spike_train::DataFrame, 
    n_neurons::Int64;
    delta::Int64=5, 
    n_bins::Int64=1,        
    duration::Float64=2000.0)
    
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

#### GC Raster plots only ####
gcs = Dict(
    "DG" => [0, 500],
)

for run_ ∈ 1:n_runs
    for i ∈ 1:length(labels)
        spikes = load_spike_files(patterns, labels[i]*"-$run_", populations)

        # CREATE RASTER PLOTS
        for p ∈ unique(spikes.Pattern)
            stimin = spikes[(spikes.Population .== "PP") .& (spikes.Pattern .== p), :]
            plots = []
            append!(plots, [raster_plot(stimin; ylab="PP")])
            for pop ∈ keys(gcs)
                lb, ub = populations[pop]
                popspikes = spikes[(spikes.Population .== pop) .& (spikes.Pattern .== p),:]
                append!(plots, [raster_plot(popspikes; xlab="", ylab=pop)])
            end
            fig = plot(reverse(plots)..., layout=grid(2, 1, heights=[0.8, 0.2]), size=(500, 400), dpi = 300)
            savefig(fig, "figures/raster-plots/raster-gconly--"*string(p)*"-"*labels[i]*"-$run_"*".png")
        end
    end 
end

#-------PLOT LFPS
l = @layout [grid(2,1, heights=[0.8, 0.2])]
theta_lfp_pp = plot(lfp_save_pp["theta"][1], 
            color = :black,
            xlabel = "Time (ms)",
            ylabel = "PP",
            label = nothing,
)          

theta_lfp_gc = plot(lfp_save_gc["theta"][1], 
            color = :black,
            ylabel = "DG",
            label = nothing,
)

sumfig = plot(theta_lfp_gc, theta_lfp_pp, layout=l, size=(500, 400), dpi=300)
savefig(sumfig, "figures/lfp/lfps_gcs_pp.png")

#-------POWER SPECTRA 
restruct_gclfps = zeros((duration, length(labels)))
restruct_pplfps = zeros((duration, length(labels)))

for i ∈ 1:length(labels)
    restruct_gclfps[:,i] = lfp_save_gc[labels[i]][n_runs]
    restruct_pplfps[:,i] = lfp_save_pp[labels[i]][n_runs]
end

# Calculate power spectra
S = spectra(restruct_gclfps, 500, duration; tapering=hamming) #TODO: learn what this 'tapering' is.

# Calculate power spectra
S = spectra(restruct_gclfps, 500, duration; tapering=hamming) #TODO: learn what this 'tapering' is.

# SUM LOW GAMMA POWER FOR 30Hz-60Hz
norm_gamma_pwrL = OrderedDict("theta"=>0.0, "ftheta"=>0.0, "alpha"=>0.0, "beta"=>0.0, "gamma"=>0.0)
for freq ∈ 1:length(labels)
    norm_gamma_pwrL[labels[freq]] =sum(S.y[30:39,freq])
end

# SUM MED GAMMA POWER FOR 30Hz-60Hz
norm_gamma_pwrM = OrderedDict("theta"=>0.0, "ftheta"=>0.0, "alpha"=>0.0, "beta"=>0.0, "gamma"=>0.0)
for freq ∈ 1:length(labels)
    norm_gamma_pwrM[labels[freq]] =sum(S.y[40:59,freq])
end

# SUM HIGH GAMMA POWER FOR 60-200Hz
norm_gamma_pwrH = OrderedDict("theta"=>0.0, "ftheta"=>0.0, "alpha"=>0.0, "beta"=>0.0, "gamma"=>0.0)
for freq ∈ 1:length(labels)
    norm_gamma_pwrH[labels[freq]] =sum(S.y[60:200,freq])
end

CSV.write("figures/lfp/gamma_pwr_low.csv", norm_gamma_pwrL)
CSV.write("figures/lfp/gamma_pwr_med.csv", norm_gamma_pwrM)
CSV.write("figures/lfp/gamma_pwr_high.csv", norm_gamma_pwrH)

# PLOT GAMMA POWER BAR GRAPH 
gamma_pwr = bar(vec(freqlabels), [[norm_gamma_pwrM[i] for i ∈ labels] [norm_gamma_pwrL[i] for i ∈ labels] [norm_gamma_pwrH[i] for i ∈ labels]], 
                        labels = ["Low "*L"\gamma" "Med. "*L"\gamma" "High "*L"\gamma"], 
                        color = [:grey75 :grey50 :black], xlabel = "Frequency Band", ylabel = L"\gamma"*" Power", 
                        xtickfont=font(12), size = (350,350), legend = :topleft, dpi=300)
savefig(gamma_pwr, "figures/lfp/gamma_pwr_bar.png")


# PLOT GAMMA POWER LINE GRAPH
gamma_pwr_line = plot(vec(freqlabels), [[norm_gamma_pwrL[i] for i ∈ labels] [norm_gamma_pwrM[i] for i ∈ labels] [norm_gamma_pwrH[i] for i ∈ labels]], 
                        labels = ["Low "*L"\gamma" "Med. "*L"\gamma" "High "*L"\gamma"], xtickfont=font(12),
                        color = [:grey75 :grey50 :black], ylims = (0, 40), linewidth = 2,  xlabel = "Input Frequency Band", ylabel = L"\gamma"*" Power", 
                        size = (350,350), legend = :topleft, dpi=300)
savefig(gamma_pwr_line, "figures/lfp/gamma_pwr_line.png")

# PLOT RAW SPECTRA
powerfig = plot(S, fmax = 50, label = freqlabels, legend = :topright, ylabel = "Power "*L"(\mu V^2)", dpi=300)
savefig(powerfig, "figures/lfp/powerspectra_gcs.png")

# PLOT HIGH FREQ RAW SPECTRA SEPERATELY
powerfigstack = plot(S, layout=grid(5,1, heights=[0.2,0.2,0.2,0.2,0.2]), xticks = (collect(0:5:50)), 
                xlims = (0,50), ylims=(0,35),label = freqlabels, linewidth = 3,
                legend = :topright, ylabel = "Power "*L"(\mu V^2)", 
                xlabel = "Frequency (Hz)", size=(500, 800), c=:black, grid=:none, dpi=300)
savefig(powerfigstack, "figures/lfp/powerspectra_gcs_stacked.png")

# PLOT FREQUENCY-TIME HEATMAP
A = TFamplitude(lfp_save_gc["theta"][n_runs], 400, duration; fmax=40)
freqtime = heatmap(A.y, xlabel = "Time (ms)", ylabel = "Frequency (Hz)", dpi=300)
savefig(freqtime, "figures/lfp/gcs_freqtime.png")



