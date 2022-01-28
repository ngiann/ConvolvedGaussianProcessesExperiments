# 27 Jan 2022
#
# Commands executed for getting results for synthetic data.
#
# This is a simulation with synthetic data.
# The purpose of the simulation is to see whether it is possible to recover
# the mass for a set of light curves that do not overlap.
#
# Here the lightcurves are even wider apart. On average there is 1 observation per unit time.
# There 4 lightcurves.
#
# Run on hephaistos
#
# ~/julia-1.7.1/bin/julia -O3 -p64
#

@everywhere using ConvolvedGaussianProcesses, ProgressMeter, Suppressor

using Printf, MiscUtil, ADDatasets, TransferFunctions, JLD2, UnicodePlots, PyPlot, Random, Distributions, ConvolvedKernel, LinearAlgebra

include("simulatenonoverlapping.jl")

tobs, yobs, σobs, tfarray, lambda = simulatenonoverlapping(N=50, mass=2e8, σ=1.0, seed=2);

PyPlot.savefig("lightcurves.svg")
PyPlot.savefig("lightcurves.png")


masses = collect(logrange(1e5, 1e10, 128));

T = [PhysicalTransferFunctions(mass = m, eddingtonfraction = 10.0, wavelengths = lambda) for m in masses]


# warmup

for i in 1:2

    outparallel = @showprogress pmap(tfarray->(@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname="matern32", tfarray=tfarray, iterations=2, numberofrestarts=1, ρmax=1.0)), T)

end

# proper run

outparallel = @showprogress pmap(tfarray->(@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname="matern32", tfarray=tfarray, iterations=2500, numberofrestarts=2, ρmax=10.0)), T)

# save results

prob = getprobabilities(outparallel)

JLD2.save("experiment6.jld2", "out", outparallel, "prob", prob, "lambda", lambda, "masses", masses, "eddingtonfraction", 10.0, "truemass", 2e8)

UnicodePlots.barplot(masses, prob)

figure()
plot(masses, prob, "o-"); xscale(:log)
plot(2e8*ones(5), collect(LinRange(minimum(prob),maximum(prob), 5)), "r--")
PyPlot.title("synthetic, non-overlapping time series, true mass is 2e8 at dashed line")
PyPlot.savefig("posteriormass.svg")
PyPlot.savefig("posteriormass.png")


# plot best fit

figure()

bestmass = masses[argmax(prob)]

tfarray = PhysicalTransferFunctions(mass=bestmass, eddingtonfraction=10.0, wavelengths=lambda);

fmin, pred = convolvedgp(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname="matern32", iterations=10000, ρmax=10.0,  tfarray = tfarray);

xtest = collect(LinRange(0.0, 360.0, 2000))

μpred, σpred = pred(xtest)

colours = ["b", "orange", "g", "m"]

for i in 1:length(tobs)
    plot(tobs[i], yobs[i], "o", color=colours[i], markeredgecolor="k")
end

for i in 1:length(tobs)
    fill_between(xtest, μpred[i] .+ σpred[i], μpred[i] .- σpred[i], color=colours[i],  alpha=0.1)
end

for i in 1:length(tobs)
    plot(xtest, μpred[i], color=colours[i])
end

PyPlot.title(@sprintf("best fit for mass=%e", bestmass))

PyPlot.savefig("bestfit.svg")
PyPlot.savefig("bestfit.png")
