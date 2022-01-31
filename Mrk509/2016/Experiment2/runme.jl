
# load data
lambda, tobs, yobs, σobs = readdataset(; source = "Mrk509_2017")

# specify kernels to try
K = ["matern12", "matern32", "rbf"]

# specify eddington fractions to try
E = [10.0, 20.0, 30.0]

# make all combinations of K and E
comb = vec(collect(Base.Iterators.product(K,E)))

# run all experiments
map(X->runexperiment(lambda = lambda[1:3], tobs = tobs[1:3], yobs = yobs[1:3], σobs = σobs[1:3], objectname = "Mrk509_2016", kernelname =  X[1], ef = X[2]), comb);
