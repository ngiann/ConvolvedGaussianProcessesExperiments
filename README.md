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

All experiments with real data use the same kernel "matern32".
As a reminder: for each object we try out a set of mass-eddington fraction combinations.
For each combination, we perform a 5-fold cross-validation and use the result to calculate a fitness value (i.e. log-likeliohood) that tells us how well the combination in question fits the data.
Once we have collected the fitnesses of all combinations, we calculate the posterior probability of the combinations. The posterior tells us how likely each combination is relatively to all other combinations. 

As an example, we look at the results of Mrk279_2017. 
We change in the corresponding folder, start Julia and load the data:
```
using JLD2, ADDatasets, PyPlot, MiscUtil, TransferFunctions
lambda, = readdataset(source = "Mrk279_2017"); # load only wavelengths
```

The folder contains a file named after the object source and extension "jld2". For the case of Mrk279_2017 the file is called "Mrk279_2017.jld2". We load this file:
```
@load "Mrk279_2017.jld2" # warnings may appear
```

Loading the file introduces a number of variables in the workspace that will be listed in the REPL once the command has been executed. The variables that we will be looking at are `masses`, `accretions`, `centroids` and `posterior`. These are all matrices with dimensions *(number of candidate masses)×(number of candidate eddington fractions)*.

Matrix `posterior` contains the posterior probability of each combination of mass and eddingtonfraction. There are 64 candidate masses (64 rows) and 3 candidate eddington fractions (3 columns) To inspect the posterior of e.g. 5th candidate mass and 2nd candidate eddington fraction, do: 
```
# create matrix of eddington fractions because they are not saved in jld2 file.
edfractions = kron([1;5;10]', ones(64)) 

i,j = 5, 2 # specify combination

posterior[i,j] # returns posterior probability
masses[i,j] # returns the mass of combination
accretions[i,j] # returns accretion rate of combination
centroids[i,j] # returns centroids per wavelength in same order as lambda when loading dataset
edfractions[i,j] # return eddington fraction of combination
```

To find the most likely combination and look at its mass, accretion and centroids, we do:
```
bestindex = argmax(posterior)
masses[bestindex]
accretions[bestindex]
centroids[bestindex]
edfractions[bestindex]
```

To plot the marginal posterior for mass, we do:
```
figure()
plot(masses[:,1], sum(posterior, dims=2), "o-", label="first"); xscale("log")
```

To plot the mass posterior for the j-th eddington fraction, we do:
```
figure()
plot(masses[:,1], posterior[:,j], "o-"); xscale("log")
```

To plot the delay posterior at k-th wavelength, for the j-th eddington fraction, we do:
```
k=1
figure()
plot([c[k] for c in centroids[:,j]], posterior[:,j], "o-"); xscale("log")
``` 

To plot the fit for the most likely combination, we load:
```
@load "predictions.jld2"
```

This introduces that variables `xtest`, `μ` and `σ` in the workspace. `xtest` holds the prediction times. `μ` and `σ` are each an array of arrays. The outer dimension goes over the wavelengths (as specified in `lambda`) and the inner array holds the prediction at the given times `xtest`.

Let us plot the predictions over the data:
```
plotdataset(source = "Mrk279_2017");
clrs = ["blue", "orange", "green", "red"]

for (index, λ) in enumerate(lambda)
  plot(xtest, μ[index], color=clrs[index])
  fill_between(xtest, μ[index] .- σ[index], μ[index] .+ σ[index], color=clrs[index], alpha=0.3)
end

```
