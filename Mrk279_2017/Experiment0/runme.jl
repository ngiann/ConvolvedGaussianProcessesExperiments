#-------------------------------------------------------------------------------
# SPECIFY DATA
#-------------------------------------------------------------------------------

SOURCE = "Synth"


lambda, tobs, yobs, σobs, tfarray = simulatedatafromgp(N=50, Tmax=100, mass=1e8, accretion=1.0, kernelname = "matern32", seed=1, σ=0.2)

lambda = lambda[1:3]
tobs   = tobs[1:3]
yobs   = yobs[1:3]
σobs   = σobs[1:3]


#-------------------------------------------------------------------------------
# SPECIFY TRANSFER FUNCTIONS
#-------------------------------------------------------------------------------

# specify kernel
kernelname = "matern32"

# specify physical parameters
masses     = collect(logrange(1e6, 1e10, 30))
accretions = collect(logrange(0.01, 3,   30))

# create combinations of transfer functions
TF = [PhysicalTransferFunctions(mass=m, accretion=a, wavelengths=lambda) for m in masses, a in accretions]


#-------------------------------------------------------------------------------
# RUN EXPERIMENTS
#-------------------------------------------------------------------------------

# warmup
runexperiment("warmup_deleteme"; transferFunctions = TF[1:min(nworkers(), length(TF))], tobs = tobs, yobs = yobs, σobs = σobs, kernelname = kernelname, fs = 100, ρmax=1.0, iterations=2)

# run all experiments

EXPERIMENT = SOURCE

savedfile = runexperiment(EXPERIMENT; transferFunctions = TF, tobs = tobs, yobs = yobs, σobs = σobs, kernelname = kernelname, fs = 100, ρmax = 9.0)


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

end
