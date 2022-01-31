# Get results for Mrk509_2016

# load data

SOURCE = "Mrk509_2016"

TF = ShakuraBoxTransferFunctions

lambda, tobs, yobs, σobs = readdataset(; source = SOURCE)

# specify kernels to try
K = ["matern12", "matern32", "rbf"]

# specify eddington fractions to try
E = [10.0, 20.0, 30.0]

# make all combinations of K and E
comb = vec(collect(Base.Iterators.product(K, E)))

# run all experiments
map(X -> runexperiment(transferFunctions = TF, lambda = lambda, tobs = tobs, yobs = yobs, σobs = σobs, objectname = SOURCE, kernelname =  X[1], ef = X[2]), comb);
