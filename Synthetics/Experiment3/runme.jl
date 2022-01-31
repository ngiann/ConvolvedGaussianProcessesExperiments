# Non-overlapping lighcurves experiments
#
# This is a simulation with synthetic data.
# The purpose of the simulation is to see whether it is possible to recover
# the mass for a set of light curves that do not overlap.
# A special script for generating the data has been written called "simulatenonoverlapping.jl".

SOURCE = "NON_OVERLAPPING"

TF = PhysicalTransferFunctions

include("simulatenonoverlapping.jl")

tobs, yobs, σobs, _tfarray_, lambda = simulatenonoverlapping(N=25, mass=2e8, σ=1.0,seed=2);

# specify kernels to try
K = ["matern12", "matern32", "rbf"]

# specify eddington fractions to try
E = [10.0, 20.0, 30.0]

# make all combinations of K and E
comb = vec(collect(Base.Iterators.product(K, E)))

# run all experiments
map(X -> runexperiment(transferFunctions = TF, lambda = lambda, tobs = tobs, yobs = yobs, σobs = σobs, objectname = SOURCE, kernelname =  X[1], ef = X[2]), comb);
