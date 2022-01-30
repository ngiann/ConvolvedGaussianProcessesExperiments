@everywhere using ConvolvedGaussianProcesses, ProgressMeter, Suppressor
using Printf, MiscUtil, ADDatasets, TransferFunctions, JLD2


function runexperiment(; lambda = lambda, tobs = tobs, yobs = yobs, σobs = σobs, objectname = objectname)


    colourprint(@sprintf("Started working on object |%s|\n", objectname))


    ###########################
    # set up parameter ranges #
    ###########################

    masses = collect(logrange(1e5, 1e10, 64))


    ##############
    # warmup run #
    ##############

    for _ in 1:3

        # create candidate transfer functions
        Φ = [PhysicalTransferFunctions(mass = m, eddingtonfraction = 10.0, wavelengths = lambda) for m in masses]

        out = @showprogress pmap(tfarray-> (@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname="matern32", tfarray=tfarray, iterations=10, numberofrestarts=1, ρmax=1.0)), Φ[1:nworkers()])

        JLD2.save("deleteme.jld2", "masses", masses, "eddingtonfraction", 10.0, "out", out)

    end

    run(`rm deleteme.jld2`)


    ##############
    # proper run #
    ##############

    for kernelname in ["matern12", "matern32", "rbf"]

        for EF in [10.0, 20.0, 30.0]

            colourprint(@sprintf("Running %s with kernel=%s and eddingtonfraction=%d\n", objectname, kernelname, EF), foreground=:red, bold=true)

            # create candidate transfer functions
            Φ = [PhysicalTransferFunctions(mass = m, eddingtonfraction = EF, wavelengths = lambda) for m in masses]

            out = @showprogress pmap(tfarray->(@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname=kernelname, tfarray=tfarray, iterations=3500, numberofrestarts=1, ρmax=20.0)), Φ)

            filename = @sprintf("%s_EF_%d_%s.jld2", objectname, Int(EF), kernelname)

            JLD2.save(filename, "objectname", objectname, "masses", masses, "eddingtonfraction", EF, "out", out, "posterior", getprobabilities(out))

            @printf("Saved results in %s\n", filename)

        end

    end


    colourprint(@sprintf("Finished working on object |%s|\n", objectname))


end
