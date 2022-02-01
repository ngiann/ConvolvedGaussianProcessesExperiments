# Get results for Mrk509_2016

# load data

# Purpose of experiment is to see whether a finer grid of masses
# will produce more "meaningfull" posterior mass distribution

SOURCE = "Mrk509_2016"

TF = PhysicalTransferFunctions

lambda, tobs, yobs, σobs = readdataset(; source = SOURCE)

# specify kernels to try
K = ["matern32"]

# specify eddington fractions to try
E = [10.0, 20.0, 30.0]

# make all combinations of K and E
comb = vec(collect(Base.Iterators.product(K, E)))

# run all experiments
map(X -> runexperiment(numberofmasses = 256, transferFunctions = TF, lambda = lambda, tobs = tobs, yobs = yobs, σobs = σobs, objectname = SOURCE, kernelname =  X[1], ef = X[2]), comb);
