# Get results for Synthetic data

# load data

SOURCE = "SYNTH"

TF = PhysicalTransferFunctions

lambda, tobs, yobs, σobs, tfarray = simulatedatafromgp(N=50, Tmax=100, mass=1e8, eddingtonfraction=10, kernelname = "matern32", seed=1)

# specify kernels to try
K = ["matern12", "matern32", "rbf"]

# specify eddington fractions to try
E = [10.0, 20.0, 30.0]

# make all combinations of K and E
comb = vec(collect(Base.Iterators.product(K, E)))

# run all experiments
map(X -> runexperiment(transferFunctions = TF, lambda = lambda, tobs = tobs, yobs = yobs, σobs = σobs, objectname = SOURCE, kernelname =  X[1], ef = X[2]), comb);
