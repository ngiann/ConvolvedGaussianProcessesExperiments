
# Purpose of this script is check out a particular model-misspecification scenario:
# This scenario is the case where one of the lightcurves does not bear any connection
# to the others. We expect the model in this case to set the scaling coefficient of
# the irrelevant lightcurve to zero and still recover the true physical parameters.


#-------------------------------------------------------------------------------
# SPECIFY DATA
#-------------------------------------------------------------------------------

SOURCE = "SYNTH"

truemass      = 1e8
trueedfrac    = 20
trueaccretion = TransferFunctions.massaccretionfunction(bhm = truemass, edfrac=trueedfrac, eta=0.1)

lambda, tobs, yobs, σobs, = simulatedatafromgp(mass=truemass, accretion=trueaccretion, Tmax = 100, N=60, seed = 1)


# generate one more dataset with different seed

lambdaX, tobsX, yobsX, σobsX, = simulatedatafromgp(mass=truemass, accretion=trueaccretion, Tmax = 100, N=60, seed = 10101)

# copy the 4th lightcurve to the original dataset
tobs[4], yobs[4], σobs[4] = tobsX[4], yobsX[4], σobsX[4]

#-------------------------------------------------------------------------------
# SPECIFY TRANSFER FUNCTIONS
#-------------------------------------------------------------------------------

# specify kernel
kernelname = "matern32"

# specify physical parameters
masses     = collect(logrange(1e6, 1e10, 32))
efractions = [10.0; 20.0; 30.0]

# create combinations of transfer functions
TF = [PhysicalTransferFunctionsEddington(mass=m, eddingtonfraction=ef, wavelengths=lambda) for m in masses, ef in efractions]



#-------------------------------------------------------------------------------
# RUN EXPERIMENTS
#-------------------------------------------------------------------------------

# warmup

for _ in 1:2

    runexperiment("warmup_deleteme"; transferFunctions = TF[1:min(nworkers(), length(TF))], tobs = tobs, yobs = yobs, σobs = σobs, kernelname = kernelname, fs = 100, ρmax=1.0, iterations=2)

    run(`rm warmup_deleteme.jld2`)

end


# run all experiments

EXPERIMENT = SOURCE

savedfile = runexperiment(EXPERIMENT; transferFunctions = TF, tobs = tobs, yobs = yobs, σobs = σobs, kernelname = kernelname, fs = 100, ρmax = 20.0)


#-------------------------------------------------------------------------------
# PLOTTING
#-------------------------------------------------------------------------------

let

    data = load(savedfile)


    #---------------------#
    # plot mass posterior #
    #---------------------#

    figure()

    title(EXPERIMENT)

    subplot(211)

    title(@sprintf("true mass is %e"), truemass)

    massmarginal = vec(sum(reshape(data["posterior"], size(TF)), dims=2))

    plot(masses, massmarginal, "o-"); xscale("log")

    subplot(212)

    title(@sprintf("true edfraction is %e"), trueedfrac)

    accmarginal = vec(sum(reshape(data["posterior"], size(TF)), dims=1))

    plot(efractions, accmarginal, "o-"); xscale("log")

    savefig("posteriors.svg"); savefig("posteriors.png")


    #---------------------------------------------#
    # fit data with most likely transfer function #
    #---------------------------------------------#

    bestTF = data["besttransferfunctionarray"]

    display(bestTF)

    _, pred = convolvedgp(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname=kernelname, tfarray=bestTF, iterations=3000, numberofrestarts = 1, ρmax=20.0)

    xtest = LinRange(minimum(map(minimum, tobs)), maximum(map(maximum, tobs)), 2500)

    μ, σ = pred(xtest)


    #----------------------#
    # plot most likely fit #
    #----------------------#

    figure()

    title(EXPERIMENT)

    clr = ["b", "orange", "g", "r", "c", "y", "m", "lightpink"]

    for i in 1:length(lambda)

        plot(tobs[i], yobs[i], "o", color = clr[i], markeredgecolor="k")

        fill_between(xtest, μ[i] .+ σ[i], μ[i] .- σ[i], color=clr[i],  alpha=0.2)

        plot(xtest, μ[i], "k-", linewidth=2)

    end

    savefig("bestfit.svg"); savefig("bestfit.png")

end
