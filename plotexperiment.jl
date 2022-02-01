using MiscUtil, PyPlot, JLD2, Statistics, StatsFuns, ConvolvedGaussianProcesses, ADDatasets, TransferFunctions, Printf, LinearAlgebra, DelimitedFiles

function plotexperiment(; filenames = filenames,
                          lambda = lambda, tobs = tobs, yobs = yobs, σobs = σobs)

    #---------------------------------------#
    #          Load saved results           #
    #---------------------------------------#

    data = map(load, filenames)

    masses = data[1]["masses"]

    # TODO : make sure all data contain the same masses

    #--------------------------------------------------------#
    # Plot posterior probabilities of each model             #
    # By model we mean eddingtonfraction-kernel combination  #
    #--------------------------------------------------------#

    figure(1); cla()

    for (d, f) in zip(data, filenames)

        plot(masses, d["posterior"], "-", label = f); xscale("log")

    end

    legend()

    savefig("massposterior.svg")

    #---------------------------------------#
    # calculate posterior model probability #
    #---------------------------------------#


    numberoffolds = length(data[1]["out"][1])

    @printf("There are %d number of folds\n", numberoffolds)

    p = zeros(length(data))

    for i in 1:length(masses), j in 1:numberoffolds

        aux = vec([d["out"][i][j] for d in data])

        p += exp.(aux .- logsumexp(aux))

    end

    p = p / (length(masses) * length(data[1]["out"][1]))

    @show sum(p)



    #----------------------------------------------------------------------------
    # Plot best fit for each combination of kernel and eddington fraction
    #----------------------------------------------------------------------------

    firstfigure = true

    # create an object type axis before we enter the loop below so that we can
    # set this object inside the loop
    fig = figure(); ax = axis(); close(fig)

    for d in data

        local kernelname, ef = d["kernelname"], d["eddingtonfraction"]

        namestring = @sprintf("%s_%s_%.2f", d["objectname"], kernelname, ef)


        # find index of most likely mass
        bestmass_index = argmax(d["posterior"])

        # get most likely mass, i.e. mode of posterior mass distribution
        bestmass = d["masses"][bestmass_index]

        @printf("best mass for %s is %e\n", namestring, bestmass)

        tfarray = PhysicalTransferFunctionsEddington(mass = bestmass, eddingtonfraction = ef, wavelengths = lambda)

        _, pred = convolvedgp(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname=kernelname, tfarray=tfarray, iterations=3000, numberofrestarts = 1, ρmax=20.0);

        xtest = LinRange(minimum(map(minimum, tobs)), maximum(map(maximum, tobs)), 2500)

        μ, σ = pred(xtest)


        figure()

        title(namestring)

        clr = ["b", "orange", "g", "r", "c", "y", "k", "m"]

        for i in 1:length(lambda)

            plot(tobs[i], yobs[i], "o", color = clr[i], markeredgecolor="k")

            fill_between(xtest, μ[i] .+ σ[i], μ[i] .- σ[i], color=clr[i],  alpha=0.2)

            plot(xtest, μ[i], "k-", linewidth=2)

            # plot(xtest, μ[i], "--", linewidth=1, color = clr[i], label = @sprintf("%d", lambda[i]))

        end

        if firstfigure
            firstfigure = false
            ax = axis()
        else
            axis(ax)
        end

        savefig(namestring * "_bestfit.svg")

        # legend()

    end

    #------------------------------------------

    writedlm("modelposterior.txt", [filenames p])

    #------------------------------------------

    return p

end
