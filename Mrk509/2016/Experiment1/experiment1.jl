# 29 Jan 2022
#
# Commands executed for getting results for Mrk509_2016
# using updated physical transfer functions that take eddingtonfraction
# instead of accretion rate
#
# Run on hephaistos
#
# ~/julia-1.7.1/bin/julia -O3 -64
#

@everywhere using ConvolvedGaussianProcesses, ProgressMeter, Suppressor
using Printf, MiscUtil, ADDatasets, TransferFunctions, JLD2

# load data
lambda, tobs, yobs, σobs = readdataset(source="Mrk509_2016");

# set up parameter ranges
masses = collect(logrange(1e5, 1e10, 64));

# warmup
for _ in 1:3

    # create candidate transfer functions
    Φ = [PhysicalTransferFunctions(mass = m, eddingtonfraction = 10.0, wavelengths = lambda) for m in masses]

    outphys = @showprogress pmap(tfarray-> (@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname="matern32", tfarray=tfarray, iterations=10, numberofrestarts=1, ρmax=1.0)), Φ[1:nworkers()]);

    JLD2.save("deleteme.jld2", "masses", masses, "eddingtonfraction", 10.0, "out", outphys)

end

run(`rm deleteme.jld2`)

for kernelname in ["matern12", "matern32", "rbf"]

    for EF in [10.0, 20.0, 30.0]

        colourprint(@sprintf("Running Mrk509_2016 with kernel=%s and eddingtonfraction=%d\n", kernelname, EF),foreground=:red, bold=true)

        Φ = [PhysicalTransferFunctions(mass = m, eddingtonfraction = EF, wavelengths = lambda) for m in masses]

        # proper run and save results
        outphys = @showprogress pmap(tfarray->(@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname=kernelname, tfarray=tfarray, iterations=3500, numberofrestarts=1, ρmax=20.0)), Φ);

        filename = @sprintf("Mrk509_2016_physical_exp1_EF_%d_%s.jld2", Int(EF), kernelname)

        JLD2.save(filename, "masses", masses, "eddingtonfraction", EF, "out", outphys, "posterior", getprobabilities(outphys))

        @printf("Saved results in %s\n", filename)

    end

end
