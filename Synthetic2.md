# Synthetic 2

Purpose of this script is to investigate how the model will behave when one of the lightcurves in the dataset bears no connection to the other lightcurves, i.e. the lightcurve in question does not come from the same latent signal as the others.
Our expectation is that the model should set the scaling coefficient of the lightcurve in question to zero and still recover the correct parameters.

Code and results for experiment [here](Synthetics/Experiment2/).

## Mass and EF posterior

![Synth1_posterior_mass](Synthetics/Experiment2/posteriors.svg)


## Most likely fit

![Synth1 best_model_fit](Synthetics/Experiment2/bestfit.svg)
