# 28 Jan 2022
#
# Commands executed for getting results for synthetic data.
#
# Same as experiment1.jl but with more noise put on the data:
# check σ argument when calling "simulatedatafromgp"
#
# Run on laptop
#
# ~/julia-1.7.1/bin/julia -O3 -p16
#

@everywhere using ConvolvedGaussianProcesses, ProgressMeter, Suppressor

using Printf, MiscUtil, ADDatasets, TransferFunctions, JLD2, UnicodePlots, PyPlot


tobs, yobs, σobs, tfarray = simulatedatafromgp(N=50, Tmax=100, mass=2e8, σ=1.0);

PyPlot.savefig("lightcurves.svg")
PyPlot.savefig("lightcurves.png")

lambda = [4300.0; 5700.0; 6200.0; 7000.0]

masses = collect(logrange(1e5, 1e10, 64));

T = [PhysicalTransferFunctions(mass = m, eddingtonfraction = 10.0, wavelengths = lambda) for m in masses]


# warmup

for i in 1:2

    outparallel = @showprogress pmap(tfarray->(@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname="matern32", tfarray=tfarray, iterations=2, numberofrestarts=1, ρmax=1.0)), T)

end

# proper run

outparallel = @showprogress pmap(tfarray->(@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname="matern32", tfarray=tfarray, iterations=2500, numberofrestarts=1, ρmax=10.0)), T)

# save results

prob = getprobabilities(outparallel)

JLD2.save("experiment2.jld2", "out", outparallel, "prob", prob, "lambda", lambda, "masses", masses, "eddingtonfraction", 10.0, "N", 50, "Tmax", 100, "truemass", 2e8, "σ", 0.2)

UnicodePlots.barplot(masses, prob)

figure()
plot(masses, prob, "o-"); xscale(:log)
plot(2e8*ones(5), collect(LinRange(minimum(prob),maximum(prob), 5)), "r--")
PyPlot.title("synthetic, true mass is 2e8")
PyPlot.savefig("posteriormass.svg")
PyPlot.savefig("posteriormass.png")


# plot best fit

figure()

bestmass = masses[argmax(prob)]

tfarray = PhysicalTransferFunctions(mass=bestmass, eddingtonfraction=10.0, wavelengths=lambda);

fmin, pred = convolvedgp(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname="matern32", iterations=10000, ρmax=10.0,  tfarray = tfarray);

xtest = collect(LinRange(0.0, 101.0, 1000))

μpred, σpred = pred(xtest)

colours = ["b", "orange", "g", "r"]

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
