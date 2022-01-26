# 26 Jan 2022
#
# Commands executed for getting results for MRK_279_2019
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
lambda, tobs, yobs, σobs = readdataset(source="Mrk279_2019");

# set up parameter ranges
masses = logrange(1e5, 1e10, nworkers());

# warmup
for _ in 1:3

    # create candidate transfer functions
    Φ = [PhysicalTransferFunctions(mass = m, eddingtonfraction = 10.0, wavelengths = lambda) for m in masses]

    outphys = @showprogress pmap(tfarray-> (@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname="matern32", tfarray=tfarray, iterations=10, numberofrestarts=1, ρmax=1.0)), Φ[1:nworkers()]);

    JLD2.save("deleteme.jld2", "masses", masses, "eddingtonfraction", 10.0, "out", outphys)

end


for EF in [5, 10, 15, 20, 25, 30.0]

    Φ = [PhysicalTransferFunctions(mass = m, eddingtonfraction = EF, wavelengths = lambda) for m in masses]

    # proper run and save results
    outphys = @showprogress pmap(tfarray->(@suppress performcv(tarray=tobs, yarray=yobs, stdarray=σobs, kernelname="matern32", tfarray=tfarray, iterations=3000, numberofrestarts=1, ρmax=20.0)), Φ);

    filename = @sprintf("Mrk279_2019_physical_exp1_EF_%d_matern32.jld2", int(EF))

    JLD2.save(filename, "masses", masses, "eddingtonfraction", EF, "out", outphys, "posterior", getprobabilities(outphys))

    @printf("Saved results in %s\n", filename)

end
