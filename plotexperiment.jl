using MiscUtil, PyPlot, JLD2, Statistics, StatsFuns, ConvolvedGaussianProcesses, ADDatasets, TransferFunctions, Printf, LinearAlgebra

function plotexperiment(; kernelnames = ["matern12", "matern32", "rbf"], EFs = [10.0, 20.0, 30.0],
                          lambda = lambda, tobs = tobs, yobs = yobs, σobs = σobs, objectname = objectname)

    #---------------------------------------#
    #      Create names and filenames       #
    #---------------------------------------#

    combinations = [(kn, ef) for kn in kernelnames, ef in EFs]

    stringcombinations = ["EF_" * @sprintf("%d", ef) * "_" * kn for kn in kernelnames, ef in EFs]

    filenames = map(x -> @sprintf("%s_%s.jld2", objectname, x), stringcombinations)


    #---------------------------------------#
    #          Load saved results           #
    #---------------------------------------#

    data = map(load, filenames)


    #--------------------------------------------------------#
    # Plot posterior probabilities of each model             #
    # By model we mean eddingtonfraction-kernel combination  #
    #--------------------------------------------------------#

    figure(1); cla()

    for (d, c) in zip(data, stringcombinations)

        plot(d["masses"], d["posterior"], "-", label = c); xscale("log")

    end

    legend()

    savefig("massposterior.svg")

    #---------------------------------------#
    # calculate posterior model probability #
    #---------------------------------------#

    p = zeros(length(data))


    for i in 1:length(masses)
        for j in 1:length(data[1]["out"][1])

            aux = vec([d["out"][i][j] for d in data])

            p += exp.(aux .- logsumexp(aux))

        end
    end

    p = p / (length(masses) * length(data[1]["out"][1]))

    @show sum(p)



    #----------------------------------------------------------------------------
    # Plot best fit for each combination of kernel and eddington fraction
    #----------------------------------------------------------------------------

    for (d, c) in zip(data, combinations)

        local kernelname, ef = c[1], c[2]

        # find index of most likely mass
        bestmass_index = argmax(d["posterior"])

        # get most likely mass, i.e. mode of posterior mass distribution
        bestmass = masses[bestmass_index]

        @printf("best mass for %s is %e\n", c, bestmass)

        tfarray = PhysicalTransferFunctions(mass = bestmass, eddingtonfraction = ef, wavelengths = lambda)

        _, pred = convolvedgp(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname=kernelname, tfarray=tfarray, iterations=3000, numberofrestarts = 1, ρmax=20.0);

        xtest = LinRange(minimum(map(minimum, tobs)), maximum(map(maximum, tobs)), 2500)

        μ, σ = pred(xtest)

        figure()

        titlestr = @sprintf("%s_%s_%f", objectname, kernelname, ef)

        title(titlestr)

        clr = ["b", "orange", "g", "r", "c", "y", "k", "m"]

        for i in 1:length(lambda)

            plot(tobs[i], yobs[i], "o", color = clr[i], markeredgecolor="k")

            fill_between(xtest, μ[i] .+ 2*σ[i], μ[i] .- 2*σ[i], color=clr[i],  alpha=0.2)

            plot(xtest, μ[i], "k-", linewidth=1, label = @sprintf("%d", lambda[i]))
        end

        savefig(titlestr * "_bestfit.svg")

        legend()

    end

    #------------------------------------------

    return p, stringcombinations

end
