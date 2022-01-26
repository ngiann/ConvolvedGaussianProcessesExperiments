using MiscUtil, PyPlot, JLD2, Statistics, StatsFuns, ConvolvedGaussianProcesses, ADDatasets, TransferFunctions, Printf, LinearAlgebra

function plot_experiment1()

    # Filenames

    file_10_matern12 = "Mrk279_2019/Experiment1/Mrk279_2019_physical_exp1_EF_10_matern12.jld2"
    file_20_matern12 = "Mrk279_2019/Experiment1/Mrk279_2019_physical_exp1_EF_20_matern12.jld2"
    file_30_matern12 = "Mrk279_2019/Experiment1/Mrk279_2019_physical_exp1_EF_30_matern12.jld2"

    file_10_matern32 = "Mrk279_2019/Experiment1/Mrk279_2019_physical_exp1_EF_10_matern32.jld2"
    file_20_matern32 = "Mrk279_2019/Experiment1/Mrk279_2019_physical_exp1_EF_20_matern32.jld2"
    file_30_matern32 = "Mrk279_2019/Experiment1/Mrk279_2019_physical_exp1_EF_30_matern32.jld2"

    file_10_rbf = "Mrk279_2019/Experiment1/Mrk279_2019_physical_exp1_EF_10_rbf.jld2"
    file_20_rbf = "Mrk279_2019/Experiment1/Mrk279_2019_physical_exp1_EF_20_rbf.jld2"
    file_30_rbf = "Mrk279_2019/Experiment1/Mrk279_2019_physical_exp1_EF_30_rbf.jld2"

    # Load saved results

    data_10_matern12 = load(file_10_matern12)
    data_20_matern12 = load(file_20_matern12)
    data_30_matern12 = load(file_30_matern12)

    data_10_matern32 = load(file_10_matern32)
    data_20_matern32 = load(file_20_matern32)
    data_30_matern32 = load(file_30_matern32)


    data_10_rbf = load(file_10_rbf)
    data_20_rbf = load(file_20_rbf)
    data_30_rbf = load(file_30_rbf)


    # Plot posterior probabilities of each model
    # By model we mean eddingtonfraction-kernel combination

    figure(1); cla()

    masses = collect(logrange(1e5, 1e10, 64))

    plot(masses, data_10_matern12["posterior"], "r-", label = "matern12, ef=10%"); xscale("log")
    plot(masses, data_20_matern12["posterior"], "g-", label = "matern12, ef=20%"); xscale("log")
    plot(masses, data_30_matern12["posterior"], "b-", label = "matern12, ef=30%"); xscale("log")

    plot(masses, data_10_matern32["posterior"], "r--", label = "matern32, ef=10%"); xscale("log")
    plot(masses, data_20_matern32["posterior"], "g--", label = "matern32, ef=20%"); xscale("log")
    plot(masses, data_30_matern32["posterior"], "b--", label = "matern32, ef=30%"); xscale("log")


    plot(masses, data_10_rbf["posterior"], "r:", label = "rbf, ef=10%"); xscale("log")
    plot(masses, data_20_rbf["posterior"], "g:", label = "rbf, ef=20%"); xscale("log")
    plot(masses, data_30_rbf["posterior"], "b:", label = "rbf, ef=30%"); xscale("log")


    legend()

    # calculate posterior model probability
    p = zeros(9)
    for i in 1:length(data_10_matern12["out"])
        for j in 1:length(data_10_matern12["out"][i])

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

    hbars = barh(1:9, p, tick_label = labels)


    hbars[1].set_color("r")
    hbars[2].set_color("g")
    hbars[3].set_color("b")

    hbars[4].set_color("r")
    hbars[5].set_color("g")
    hbars[6].set_color("b")

    hbars[7].set_color("r")
    hbars[8].set_color("g")
    hbars[9].set_color("b")

    # Plot best fit for each combination of kernel and eddington fraction

    ##################
    # Matern 12, 10% #
    ##################

    bestmass_index = argmax(data_10_matern12["posterior"])
    bestmass = masses[bestmass_index]
    @printf("best mass for matern12, 10%% is %e\n", bestmass)

    lambda, tobs, yobs, σobs = readdataset(; source = "Mrk279_2019");

    tfarray = PhysicalTransferFunctions(mass = bestmass, eddingtonfraction = 10.0, wavelengths = lambda)

    _, pred = convolvedgp(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname="matern32", tfarray=tfarray, iterations=3000, numberofrestarts = 1, ρmax=10.0);

    xtest = 57866:0.1:57957.0

    μ, σ = pred(xtest)

    figure()
    title("Mrk279_2019, matern12, ef=10%")
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

    lambda, tobs, yobs, σobs = readdataset(; source = "Mrk279_2019");

    tfarray = PhysicalTransferFunctions(mass = bestmass, eddingtonfraction = 20.0, wavelengths = lambda)

    _, pred = convolvedgp(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname="rbf", tfarray=tfarray, iterations=3000, numberofrestarts = 1, ρmax=10.0);

    xtest = 57866:0.1:57957.0

    μ, Σ = pred(xtest)

    figure()
    title("Mrk279_2019, rbf, ef=20%")
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
