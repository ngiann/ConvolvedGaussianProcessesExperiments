#
# To run multiple kernels and eddingtonfractions for the same dataset, use sth linke:
#
# K = ["matern12", "matern32", "rbf"]
# E = [10.0, 20.0, 30.0]
# comb = vec(collect(Base.Iterators.product(K,E)))
# map(X->runexperiment(lambda = lambda, tobs = tobs, yobs = yobs, σobs = σobs, objectname = "Mrk279_2017", kernelname =  X[1], ef = X[2]), comb);
#

@everywhere using ConvolvedGaussianProcesses, ProgressMeter, Suppressor
using Printf, MiscUtil, ADDatasets, TransferFunctions, JLD2


function runexperiment(; lambda = lambda, tobs = tobs, yobs = yobs, σobs = σobs, objectname = objectname, kernelname = kernelname, ef = ef)


    colourprint(@sprintf("Started working on object |%s|\n", objectname), bold = true)

    colourprint(@sprintf("kernelname is |%s|\n", kernelname))

    colourprint(@sprintf("ef is |%f|\n", ef))

    ###########################
    # set up parameter ranges #
    ###########################

    masses = collect(logrange(1e5, 1e10, 64))


    ##############
    # warmup run #
    ##############

    for _ in 1:2

        # create candidate transfer functions
        Φ = [PhysicalTransferFunctions(mass = m, eddingtonfraction = ef, wavelengths = lambda) for m in masses]

        out = @showprogress pmap(tfarray-> (@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname=kernelname, tfarray=tfarray, iterations=10, numberofrestarts=1, ρmax=1.0)), Φ[1:nworkers()])

        JLD2.save("deleteme.jld2", "masses", masses, "eddingtonfraction", 10.0, "out", out, "kernelname", kernelname)

    end

    run(`rm deleteme.jld2`)


    ##############
    # proper run #
    ##############


    colourprint(@sprintf("Running %s with kernel=%s and eddingtonfraction=%d\n", objectname, kernelname, ef), foreground=:red, bold=true)

    # create candidate transfer functions
    Φ = [PhysicalTransferFunctions(mass = m, eddingtonfraction = ef, wavelengths = lambda) for m in masses]

    out = @showprogress pmap(tfarray->(@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname=kernelname, tfarray=tfarray, iterations=3500, numberofrestarts=1, ρmax=20.0)), Φ)

    filename = @sprintf("%s_EF_%d_%s.jld2", objectname, Int(ef), kernelname)

    JLD2.save(filename, "objectname", objectname, "masses", masses, "eddingtonfraction", ef, "out", out, "posterior", getprobabilities(out), "kernelname", kernelname)

    @printf("Saved results in %s\n", filename)



    colourprint(@sprintf("Finished working on object |%s|\n", objectname))


end
