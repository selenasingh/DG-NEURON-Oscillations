using Plots, DataFrames, CSV, LaTeXStrings

df = CSV.read("ficurve.txt", delim="\t", header=0, DataFrame)
k_unique = unique(df[:,4])

df

CSV.write("test.dat", df[:,1:2], delim="\t", header=false)

labels = ["HC", "LR", "NR"]
ls = [:solid, :dash, :dot]
p = plot(xlabel=L"pA", ylabel=L"Hz", legend=:topleft)
for i ∈ 1:3
    p = plot!(1000 .* df[df[:,4] .== k_unique[i],1], df[df[:,4] .== k_unique[i], 2], 
                c=:black, ls=ls[i], label=labels[i])
end
p

V = DataFrame()
for i ∈ 0:2
    v0 = CSV.read("vm-$i.000000.txt", delim="\t", header=0, DataFrame)
    v0[:,"ID"] .= i
    V = [V; v0]
end

plot(v0[:,1], v0[:,22])

p = plot()
for i ∈ 0:2
    p = plot!(V[V.ID .== i,1], V[V.ID .== i, 22], label="$i")
end
p
