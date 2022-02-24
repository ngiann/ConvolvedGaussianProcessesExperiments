#-------------------------------------------------------------------------------
# SPECIFY DATA
#-------------------------------------------------------------------------------

SOURCE = "Mrk509_2016"


# ❗ Same as experiment 1 but only two bands, the 1st and 3rd❗

lambda, tobs, yobs, σobs = readdataset(source = SOURCE)


lambda, tobs, yobs, σobs = lambda[[1;3]], tobs[[1;3]], yobs[[1;3]], σobs[[1;3]]

#-------------------------------------------------------------------------------
# SPECIFY TRANSFER FUNCTIONS
#-------------------------------------------------------------------------------

# specify kernel
kernelname = "matern32"

# specify physical parameters
masses     = collect(logrange(1e5, 1e10, 64))
efractions = [1.0; 5.0; 10.0]

# create combinations of transfer functions
TF = pmap(((m,ef),) -> PhysicalTransferFunctionsEddington(mass = m, eddingtonfraction = ef, wavelengths = lambda), Iterators.product(masses, efractions))

include("../../../../runmecommonpart.jl")
