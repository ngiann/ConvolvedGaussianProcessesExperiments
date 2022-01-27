# 27 Jan 2022
#
# Commands executed for getting results for synthetic data.
#
# This is a simulation with synthetic data.
# The purpose of the simulation is to see whether it is possible to recover
# the mass for a set of light curves that do not overlap.
#
# same as experiment3 but more observations per lightcurve, note argument N=50.
#
# Run on hephaistos
#
# ~/julia-1.7.1/bin/julia -O3 -64
#

@everywhere using ConvolvedGaussianProcesses, ProgressMeter, Suppressor

using Printf, MiscUtil, ADDatasets, TransferFunctions, JLD2, UnicodePlots, PyPlot, Random, Distributions, ConvolvedKernel, LinearAlgebra

include("simulatenonoverlapping.jl")

tobs, yobs, σobs, tfarray, lambda = simulatenonoverlapping(N=50, mass=2e8, σ=1.0,seed=2);

PyPlot.savefig("lightcurves.svg")
PyPlot.savefig("lightcurves.png")


masses = collect(logrange(1e5, 1e10, 128));

T = [PhysicalTransferFunctions(mass = m, eddingtonfraction = 10.0, wavelengths = lambda) for m in masses]


# warmup

for i in 1:2

    outparallel = @showprogress pmap(tfarray->(@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname="matern32", tfarray=tfarray, iterations=2, numberofrestarts=1, ρmax=1.0)), T)

end

# proper run

outparallel = @showprogress pmap(tfarray->(@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname="matern32", tfarray=tfarray, iterations=2500, numberofrestarts=1, ρmax=10.0)), T)

# save results

prob = getprobabilities(outparallel)

JLD2.save("experiment4.jld2", "out", outparallel, "prob", prob, "lambda", lambda, "masses", masses, "eddingtonfraction", 10.0, "truemass", 2e8)

UnicodePlots.barplot(masses, prob)

figure()
plot(masses, prob, "o-"); xscale(:log)
PyPlot.title("synthetic, non-overlapping time series, true mass is 2e8")
PyPlot.savefig("posteriormass.svg")
PyPlot.savefig("posteriormass.png")
