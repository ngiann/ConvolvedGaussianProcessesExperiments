# Get results for Mrk509_2016

SOURCE = "Mrk509_2016"

lambda, tobs, yobs, σobs = readdataset(; source = SOURCE)


#-------------------------------------------------------------------------------
# SPECIFY TRANSFER FUNCTIONS
#-------------------------------------------------------------------------------

# specify kernel
kernelname = "matern32"

# specify physical parameters
masses     = collect(logrange(1e6, 1e10, 64))

# create combinations of transfer functions
TF = [PhysicalTransferFunctionsEddington(mass=m, eddingtonfraction=1.0, wavelengths=lambda) for m in masses]



#-------------------------------------------------------------------------------
# RUN EXPERIMENTS
#-------------------------------------------------------------------------------

# warmup

for _ in 1:2

    runexperiment("warmup_deleteme"; transferFunctions = TF[1:min(nworkers(), length(TF))], tobs = tobs, yobs = yobs, σobs = σobs, kernelname = kernelname, fs = 100, ρmax=1.0, iterations=2)

end


# run all experiments

EXPERIMENT = SOURCE * "PhysicalTF"

savedfile = runexperiment(EXPERIMENT; transferFunctions = TF, tobs = tobs, yobs = yobs, σobs = σobs, kernelname = kernelname, fs = 200, ρmax = 20.0)


#-------------------------------------------------------------------------------
# PLOTTING
#-------------------------------------------------------------------------------

let

    data = load(savedfile)


    # plot mass posterior

    figure()

    title(EXPERIMENT)

    subplot(211)

    massmarginal = vec(sum(reshape(data["posterior"], size(TF)), dims=2))

    plot(masses, massmarginal, "o-"); xscale("log")

    subplot(212)

    accmarginal = vec(sum(reshape(data["posterior"], size(TF)), dims=1))

    plot(accretions, accmarginal, "o-"); xscale("log")

    savefig("posteriors.svg")


    # fit data with most likely transfer function

    bestTF = data["besttransferfunctionarray"]

    display(bestTF)

    _, pred = convolvedgp(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname=kernelname, tfarray=bestTF, iterations=3000, numberofrestarts = 1, ρmax=20.0)

    xtest = LinRange(minimum(map(minimum, tobs)), maximum(map(maximum, tobs)), 2500)

    μ, σ = pred(xtest)


    # plot most likely fit

    figure()

    title(EXPERIMENT)

    clr = ["b", "orange", "g", "r", "c", "y", "m", "lightpink"]

    for i in 1:length(lambda)

        plot(tobs[i], yobs[i], "o", color = clr[i], markeredgecolor="k")

        fill_between(xtest, μ[i] .+ σ[i], μ[i] .- σ[i], color=clr[i],  alpha=0.2)

        plot(xtest, μ[i], "k-", linewidth=2)

    end

    savefig("bestfit.svg")

end
