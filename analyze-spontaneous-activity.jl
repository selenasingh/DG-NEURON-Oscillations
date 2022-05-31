include("utilities.jl");
default(show=true)

duration = 200
scaling_factor = 1000/duration
patterns = 0:1
populations = Dict(
    "DG" => [0, 500],
    "BC" => [500,506],
    "MC" => [506, 521],
    "HIPP" => [521, 527]
)


spikes = load_spike_files(patterns, "HC"*"-1", populations)

dgspikes = spikes[(spikes.Population .== "DG") .& (spikes.Pattern .== 0), :]
dgids = dgspikes.Neuron |> unique

[size(dgspikes[dgspikes.Neuron .== i,:], 1) for i ∈ dgids ] |> unique

freq = mean([scaling_factor * size(dgspikes[dgspikes.Neuron .== i,:], 1) / duration for i ∈ dgids ])
scaling_factor * size(dgspikes, 1) / duration / 500

sastim = CSV.read("data/dgnetwork/NoiseStimIn-0-6-HC-1.txt", delim="\t", header=0, DataFrame)

scatter(sastim[:,1], sastim[:,2], markershape=:vline)

[size(sastim[sastim[:,2] .== i,:], 1) for i ∈ unique(sastim[:,2]) ]