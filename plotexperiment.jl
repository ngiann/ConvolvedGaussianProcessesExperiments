using MiscUtil, PyPlot, JLD2, Statistics, StatsFuns, ConvolvedGaussianProcesses, ADDatasets, TransferFunctions, Printf, LinearAlgebra

function plotexperiment(objectname)


    #---------------------------------------#
    #      Create names and filenames       #
    #---------------------------------------#

    kernelnames = ["matern12", "matern32", "rbf"]

    EFs = map(x -> @sprintf("%d", x), [10, 20, 30])

    combinations = [kn * "_" * ef for kn in kernelnames, ef in EFs]

    filenames = map(x -> @sprintf("%s_EF_%s.jld2", objectname, x), combinations)


    #---------------------------------------#
    #          Load saved results           #
    #---------------------------------------#

    data = map(load, filenames)


    #--------------------------------------------------------#
    # Plot posterior probabilities of each model             #
    # By model we mean eddingtonfraction-kernel combination  #
    #--------------------------------------------------------#

    figure(1); cla()

    for (d, c) in (data, combinations)

        plot(d["masses"], d["posterior"], "r-", label = c); xscale("log")

    end

    legend()

    return
    #---------------------------------------#
    # calculate posterior model probability #
    #---------------------------------------#

    p = zeros(length(data))

    for i in 1:length(data[1]["out"])
        for j in 1:length(data[1]["out"][i])

            aux = [data_10_matern12["out"][i][j];
                   data_20_matern12["out"][i][j];
                   data_30_matern12["out"][i][j];
                   data_10_matern32["out"][i][j];
                   data_20_matern32["out"][i][j];
                   data_30_matern32["out"][i][j];
                   data_10_rbf["out"][i][j];
                   data_20_rbf["out"][i][j];
                   data_30_rbf["out"][i][j]]

            p += exp.(aux .- logsumexp(aux))

        end
    end

    p = p / (length(data_10_matern12["out"]) * length(data_10_matern12["out"][1]))

    @show sum(p)

    # Create horizontal bar plot for joint posterior probabilities
    # of kernel and eddington fraction

    ax = figure(2); cla()

    labels = ["matern12,10%","matern12,20%","matern12,30%",
              "matern32,10%","matern32,20%","matern32,30%",
              "rbf,10%","rbf,20%","rbf,30%"]

    hbars = barh(1:length(data), p, tick_label = combinations)

    # hbars[1].set_color("r")
    # hbars[2].set_color("g")
    # hbars[3].set_color("b")
    #
    # hbars[4].set_color("r")
    # hbars[5].set_color("g")
    # hbars[6].set_color("b")
    #
    # hbars[7].set_color("r")
    # hbars[8].set_color("g")
    # hbars[9].set_color("b")

    # Plot best fit for each combination of kernel and eddington fraction

    ##################
    # Matern 32, 10% #
    ##################

    bestmass_index = argmax(data_10_matern32["posterior"])
    bestmass = masses[bestmass_index]
    @printf("best mass for matern32, 10%% is %e\n", bestmass)

    lambda, tobs, yobs, σobs = readdataset(; source = "Mrk509_2016");

    tfarray = PhysicalTransferFunctions(mass = bestmass, eddingtonfraction = 10.0, wavelengths = lambda)

    _, pred = convolvedgp(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname="matern32", tfarray=tfarray, iterations=3000, numberofrestarts = 1, ρmax=10.0);

    xtest = 57525:0.5:57700

    μ, σ = pred(xtest)

    figure()
    title("Mrk509_2016, matern32, ef=10%")
    clr = ["m", "b", "g", "r"]
    for i in 1:4
        plot(tobs[i], yobs[i], clr[i]*"o", markeredgecolor="k")

        fill_between(xtest, μ[i] .+ σ[i], μ[i] .- σ[i], color=clr[i],  alpha=0.1)

        plot(xtest, μ[i], "k-", linewidth=2, label = @sprintf("%d", lambda[i]))
    end
    # legend()

    ax = axis()

    ##################
    # RBF, 20%       #
    ##################

    bestmass_index = argmax(data_20_rbf["posterior"])
    bestmass = masses[bestmass_index]
    @printf("best mass for rbf, 20%% is %e\n", bestmass)

    lambda, tobs, yobs, σobs = readdataset(; source = "Mrk509_2016");

    tfarray = PhysicalTransferFunctions(mass = bestmass, eddingtonfraction = 20.0, wavelengths = lambda)

    _, pred = convolvedgp(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname="rbf", tfarray=tfarray, iterations=3000, numberofrestarts = 1, ρmax=10.0);

    xtest = 57525:0.5:57700

    μ, σ = pred(xtest)

    figure()
    title("Mrk509_2016, rbf, ef=20%")
    clr = ["m", "b", "g", "r"]
    for i in 1:4
        plot(tobs[i], yobs[i], clr[i]*"o", markeredgecolor="k")

        fill_between(xtest, μ[i] .+ σ[i], μ[i] .- σ[i], color=clr[i],  alpha=0.1)

        plot(xtest, μ[i], "k-", linewidth=2, label = @sprintf("%d", lambda[i]))
    end
    # legend()

    axis(ax)

    #------------------------------------------

    return p
end
