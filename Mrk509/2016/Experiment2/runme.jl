# Get results for Mrk509_2017
# Use only three first lighcurves
# In Experiment1 we saw that the lightcurve at the last (highest) wvelength could not
# be modelled and the scale parameter was set to 0.0.
# We speculate that the last lightcurve cannot be modelled as the convolution with
# the transfer function, hence our model sets its scale to 0.0 to get rid of it.

# load data
lambda, tobs, yobs, σobs = readdataset(; source = "Mrk509_2017")

# specify kernels to try
K = ["matern12", "matern32", "rbf"]

# specify eddington fractions to try
E = [10.0, 20.0, 30.0]

# make all combinations of K and E
comb = vec(collect(Base.Iterators.product(K,E)))

# run all experiments
map(X -> runexperiment(lambda = lambda[1:3], tobs = tobs[1:3], yobs = yobs[1:3], σobs = σobs[1:3], objectname = "Mrk509_2016", kernelname =  X[1], ef = X[2]), comb);
