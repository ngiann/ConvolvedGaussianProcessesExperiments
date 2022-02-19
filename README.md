# ConvolvedGaussianProcessesExperiments

Refreshed experiments

## Results
- Real data
  - [Mrk279_2017](Mrk279_2017.md)
  - [Mrk509_2016](Mrk509_2016.md)
  - [Mrk509_2017](Mrk509_2017.md)
  - Ngc2617
  - [Ark120](Ark120.md)
- Synthetic
  - [Synthetic 1](Synthetic1.md)
  - Synthetic 2 - incongruous lightcurve
  - Synthetic 3 - incongruous lightcurve

## How to read the results for real data

All experiments with real data use the same grid of candidate values and the same kernel:
```
# specify kernel
kernelname = "matern32"

# specify physical parameters
masses     = collect(logrange(1e5, 1e10, 64))
efractions = [1.0; 5.0; 10.0]
```

As an example, we look at the results of Mrk279_2017. 
We change in the corresponding folder, start Julia and load the data:
```
using JLD2, ADDatasets, PyPlot, MiscUtil, TransferFunctions
lambda, tobs, yobs, σobs = readdataset(source = "Mrk279_2017");
```

The folder contains a file named after the source and extension "jld2". For the case of Mrk279_2017 the file is called "Mrk279_2017.jld2". We load this file:
```
@load "Mrk279_2017.jld2" # warnings may appear
```

Loading the file introduces a number of variables that will be listed in the REPL once the command has been executed. The variables that we will be looking at are `masses`, `accretions`, `centroids` and `posterior`. These are all matrices with dimensions *(number of candidate masses)×(number of candidate eddington fractions)*.

Matrix `posterior` containts the posterior probability of each combination of mass and eddingtonfraction. There are 64 candidate masses, hence 64 rows, and 3 candidate eddington fractions, hence 3 columns. To inspect the posterior of e.g. 5th candidate mass and 2nd candidate eddington fraction, do: 
```
i,j = 5, 2 # specify combination
posterior[i,j] # returned posterior probability
masses[i,j] # returns the mass of combination
accretions[i,j] # returns accretion rate of combination
centroids[i,j] # returns centroids per wavelength in same order as lambda when loading dataset
```

To find the most likely combination:
```
bestindex = argmax(posterior)
masses[bestindex]
accretions[bestindex]
centroids[bestindex]
```



