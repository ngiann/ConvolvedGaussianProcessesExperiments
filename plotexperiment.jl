using MiscUtil, PyPlot, JLD2, Statistics, StatsFuns, ConvolvedGaussianProcesses, ADDatasets, TransferFunctions, Printf, LinearAlgebra

function plotexperiment(; filenames = filenames,
                          lambda = lambda, tobs = tobs, yobs = yobs, σobs = σobs)

    #---------------------------------------#
    #          Load saved results           #
    #---------------------------------------#

    data = map(load, filenames)


    #--------------------------------------------------------#
    # Plot posterior probabilities of each model             #
    # By model we mean eddingtonfraction-kernel combination  #
    #--------------------------------------------------------#

    figure(1); cla()

    for (d, f) in zip(data, filenames)

        plot(d["masses"], d["posterior"], "-", label = f); xscale("log")

    end

    legend()

    savefig("massposterior.svg")

    #---------------------------------------#
    # calculate posterior model probability #
    #---------------------------------------#

    p = zeros(length(data))


    numberoffolds = length(data[1]["out"][1])

    for i in 1:length(masses, j in 1:numberoffolds

        aux = vec([d["out"][i][j] for d in data])

        p += exp.(aux .- logsumexp(aux))

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
        bestmass = d["masses"][bestmass_index]

        @printf("best mass for %s is %e\n", c, bestmass)

        tfarray = PhysicalTransferFunctions(mass = bestmass, eddingtonfraction = ef, wavelengths = lambda)

        _, pred = convolvedgp(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname=kernelname, tfarray=tfarray, iterations=3000, numberofrestarts = 1, ρmax=20.0);

        xtest = LinRange(minimum(map(minimum, tobs)), maximum(map(maximum, tobs)), 2500)

        μ, σ = pred(xtest)

        figure()

        titlestr = @sprintf("%s_%s_%f", d["objectname"], kernelname, ef)

        title(titlestr)

        clr = ["b", "orange", "g", "r", "c", "y", "k", "m"]

        for i in 1:length(lambda)

            plot(tobs[i], yobs[i], "o", color = clr[i], markeredgecolor="k")

            fill_between(xtest, μ[i] .+ σ[i], μ[i] .- σ[i], color=clr[i],  alpha=0.2)

            plot(xtest, μ[i], "k-", linewidth=1, label = @sprintf("%d", lambda[i]))
        end

        savefig(titlestr * "_bestfit.svg")

        legend()

    end

    #------------------------------------------

    return p

end
