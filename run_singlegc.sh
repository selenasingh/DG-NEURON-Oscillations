# RUN EXPERIMENTS: stimulate one GC by one oscillating PP neuron
mkdir data
mkdir data/gcpp
nrnivmodl mods
nrniv gcppstimulate.hoc