@everywhere using ConvolvedGaussianProcesses, ProgressMeter, Suppressor, TransferFunctions, MiscUtil
using Printf, ADDatasets, JLD2
using ConvolvedKernel, Random, Distributions, LinearAlgebra, PyPlot

function runexperiment(experimentname; tobs = tobs, yobs = yobs, σobs = σobs, kernelname = kernelname, transferFunctions = transferFunctions, iterations = 3500, ρmin = 0.01, ρmax = 20.0, fs = FS, T = 800)


    colourprint(@sprintf("Started experiment |%s|\n", experimentname), bold = true)

    colourprint(@sprintf("kernelname is |%s|\n", kernelname), foreground =:light_cyan)

    colourprint(@sprintf("number of candidate transfer functions is |%d|\n", length(transferFunctions)), foreground =:light_cyan)

    colourprint(@sprintf("running for a maximum |%d| number of iterations\n",iterations), foreground =:light_cyan)

    colourprint(@sprintf("minimum, maximum length scale ρmin, ρmax set to |%f|, |%f|\n", ρmin, ρmax), foreground =:light_cyan)

    colourprint(@sprintf("fs set to |%d|\n", fs), foreground =:light_cyan)

    colourprint(@sprintf("T  set to |%f|\n",  T), foreground =:light_cyan)


    out = @showprogress pmap(tfarray->(@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname=kernelname, tfarray=tfarray, iterations=iterations, numberofrestarts=1, ρmin=ρmin, ρmax=ρmax, fs = fs, T = T)), transferFunctions)


    masses     = [mass(tf[1])      for tf in transferFunctions]
    accretions = [accretion(tf[1]) for tf in transferFunctions]
    centroids  = [[centroid(tf[i]) for i in 1:length(tf)]  for tf in transferFunctions]


    filename = @sprintf("%s.jld2", experimentname)

    JLD2.save(filename, "experimentname", experimentname,
                        "out", out,
                        "posterior", reshape(getprobabilities(out), size(out)),
                        "masses",     masses,
                        "accretions", accretions,
                        "centroids",  centroids,
                        "kernelname", kernelname,
                        "tobs", tobs,
                        "yobs", yobs,
                        "σobs", σobs,
                        "ρmin", ρmin,
                        "ρmax", ρmax,
                        "fs", fs,
                        "iterations", iterations,
                        "besttransferfunctionarray", transferFunctions[argmax(getprobabilities(out))])


    colourprint(@sprintf("Saved results in %s\n", filename), foreground =:light_cyan)

    colourprint(@sprintf("Finished experiment |%s|\n\n", experimentname), bold = true)

    return filename

end
