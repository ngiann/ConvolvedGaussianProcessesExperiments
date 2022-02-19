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




