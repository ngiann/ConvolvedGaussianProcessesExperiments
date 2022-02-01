using MiscUtil, PyPlot, JLD2, Statistics, StatsFuns, ConvolvedGaussianProcesses, ADDatasets, TransferFunctions, Printf, LinearAlgebra, DelimitedFiles


function produceimages(experimentname)

    #-------------------------------------------------#
    # Get hold of all jld2 files in present directory #
    # that start with experimentname and end in .jld2 #
    #-------------------------------------------------#

    regularexpression = Regex(@sprintf("^%s.*jld2", experimentname))

    filenames = filter(x -> occursin(regularexpression, x), readdir())

    colourprint(@sprintf("Got hold of the %d files below:\n", length(filenames)), foreground =:light_cyan)

    display(filenames)


    #---------------------------------------#
    #          Load saved results           #
    #---------------------------------------#

    data = map(load, filenames)

    masses = data[1]["masses"]

    # make sure all data contain the same masses
    @assert all([d["masses"] == masses for d in data])


    #---------------------------------------------#
    # Plot posterior probabilities of each model  #
    # loaded from each file                       #
    #---------------------------------------------#

    figure(1); cla()

    for (d, f) in zip(data, filenames)

        plot(masses, d["posterior"], "-", label = f); xscale("log")

    end

    legend()

    maximisefigure()

    savefig("massposterior.svg"); close()


    #---------------------------------------#
    # calculate posterior model probability #
    #---------------------------------------#

    numberoffolds = length(data[1]["out"][1])

    # TODO: make sure all have the same folds

    @printf("There are %d number of folds\n", numberoffolds)

    p = zeros(length(data))

    for i in 1:length(masses), j in 1:numberoffolds

        aux = vec([d["out"][i][j] for d in data])

        p += exp.(aux .- logsumexp(aux))

    end

    p = p / (length(masses) * numberoffolds)

    write("posteriors.md", producemdtable(filenames, p))

    
    #--------------------------------------------------------------#
    # Plot best fit for each model loaded from the respective file #
    #--------------------------------------------------------------#

    firstfigure = true

    # create an object type axis before we enter the loop below so that we can
    # set this object inside the loop
    fig = figure(); ax = axis(); close(fig)

    for (d,f) in (data, filenames)

        namestring = replace(f, ".jld2" => "")

        # get transfer function at maximum probability

        bestTF = d["besttransferfunctionarray"]

        display(bestTF)

        kernelname = d["kernelname"]

        tobs, yobs, σobs = d["tobs"], d["yobs"], d["σobs"]

        _, pred = convolvedgp(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname=kernelname, tfarray=bestTF, iterations=3000, numberofrestarts = 1, ρmax=20.0)

        xtest = LinRange(minimum(map(minimum, tobs)), maximum(map(maximum, tobs)), 2500)

        μ, σ = pred(xtest)


        figure()

        title(namestring)

        clr = ["b", "orange", "g", "r", "c", "y", "m", "lightpink"]

        for i in 1:length(lambda)

            plot(tobs[i], yobs[i], "o", color = clr[i], markeredgecolor="k")

            fill_between(xtest, μ[i] .+ σ[i], μ[i] .- σ[i], color=clr[i],  alpha=0.2)

            plot(xtest, μ[i], "k-", linewidth=2)

        end

        if firstfigure
            firstfigure = false
            ax = axis()
        else
            axis(ax)
        end

        # legend()

        maximisefigure()

        savefig(namestring * "_bestfit.svg")

        close()

    end

end



function maximisefigure()
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()
end

function producemdtable(names, p)

    @assert(length(names) == length(p))

    output = ""
    output *= @sprintf("|name          | probability |\n")
    output *= @sprintf("|--------------|-------------|\n")

    for n in 1:length(p)
        output *= @sprintf("|%s|%.3f|\n", names[n], p[n])
    end

    return output

end
