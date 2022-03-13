# ConvolvedGaussianProcessesExperiments

Refreshed experiments

## Results
- Real data
  - [Mrk279_2017](Mrk279_2017.md)
  - [Mrk509_2016](Mrk509_2016.md)
  - [Mrk509_2017](Mrk509_2017.md)
  - [Ngc2617](Ngc2617.md)
  - [Ark120](Ark120.md)
  - [Mgc0811](Mgc0811.md)
- Synthetic
  - [Synthetic 1](Synthetic1.md)
  - [Synthetic 1 with box functions](Synthetic1Box.md)

## How to read the results for real data


##### Brief overview

All experiments on real data use the "matern32" kernel. For each object we try out a set of mass-eddington fraction combinations. In all experiments, there are 64 candidate masses and 3 candidate eddington-fractions which result in 192 possible combinations.
For each combination we perform a 5-fold cross-validation and use the result to calculate a fitness value (i.e. log-likelihood) that tells us how well the combination in question fits the data.
Once we have collected the fitness of all combinations, we calculate the posterior probability of the combinations. The posterior probability tells us how likely each combination is relatively to all other combinations. 

##### Mrk279_2017 example

As an example, we look at the results for Mrk279_2017. 
We change in the corresponding folder (in this case `RELATIVEPATH/ConvolvedGaussianProcessesExperiments/Real/Mrk279_2017/Experiment1/`), start Julia and load the data:
```
using Printf, JLD2, ADDatasets, PyPlot, MiscUtil, TransferFunctions
lambda, = readdataset(source = "Mrk279_2017"); # load only wavelengths
```

The folder contains a file named after the object source and extension "jld2". For the case of Mrk279_2017 the file is called "Mrk279_2017.jld2". We load this file:
```
@load "Mrk279_2017.jld2" # warnings may appear
```

We also create matrix of eddington fractions because they are not saved in jld2 file.
```
edfractions = kron([1;5;10]', ones(64)) 
```
Loading the file introduces a number of variables in the workspace. We will look at the variables `masses`, `accretions`, `centroids`, `posterior` that we just loaded and `edfractions` that we just created. These  matrices have dimensions *number of candidate masses*×*number of candidate eddington fractions*:
```
map(size, [edfractions, masses, accretions, centroids, posterior])
```
We see that all matrices have 64 rows (that correspond to the 64 candidate masses) and 3 rows (that correspond to the 3 candidate eddington fractions).

Each entry in matrix `posterior` contains the posterior probability of a particular combination of mass and eddingtonfraction.
The elements in matrix `posterior` correspond to the elements in the matrices `masses`, `accretions`, `centroids`,  and `edfractions`.
This means that if we want to find out about e.g. the 5th candidate mass and 2nd candidate eddington fraction, do: 
```
i,j = 5, 2 # specify combination

posterior[i,j] # returns posterior probability
masses[i,j] # returns the mass of combination
accretions[i,j] # returns accretion rate of combination
centroids[i,j] # returns centroids per wavelength in same order as lambda when loading dataset
edfractions[i,j] # return eddington fraction of combination
```

##### Most likely combination

To find the most likely combination (i.e. highest posterior probability) and look at its mass, accretion, centroids and eddington fraction, we do:
```
bestindex = argmax(posterior)
masses[bestindex]
accretions[bestindex]
centroids[bestindex]
edfractions[bestindex]
```

##### Mass posterior

To plot the marginal posterior for mass (over all Eddington fractions), we do:
```
figure()
plot(masses[:,1], sum(posterior, dims=2), "o-", label="first"); xscale("log")
```

To plot the mass posterior for the j-th eddington fraction (conditional posterior), we do:
```
figure()
plot(masses[:,1], posterior[:,j], "o-"); xscale("log")
```

##### Centroid posterior ⚠️ This needs re-thinking, don't trust for now

To plot the centroid posterior at k-th wavelength, for the j-th eddington fraction, we do:
```
k=1
figure()
plot([c[k] for c in centroids[:,j]], posterior[:,j], "o-"); xscale("log")
``` 

##### Predictions for most likely combination

To plot the fit for the most likely combination, we load:
```
@load "predictions.jld2"
```

This introduces variables `xtest`, `μ` and `σ` in the workspace. `xtest` holds the prediction times. `μ` and `σ` are each an array of arrays. The outer dimension goes over the wavelengths (as specified in `lambda`) and the inner array holds the prediction at the given times `xtest`.

Let us plot the predictions over the observed data:
```
plotdataset(source = "Mrk279_2017");
clrs = ["blue", "orange", "green", "red", "magenta", "brown"]

for (index, λ) in enumerate(lambda)
  plot(xtest, μ[index], color=clrs[index], label=@sprintf("%f", λ))
  fill_between(xtest, μ[index] .- σ[index], μ[index] .+ σ[index], color=clrs[index], alpha=0.3)
end

```
