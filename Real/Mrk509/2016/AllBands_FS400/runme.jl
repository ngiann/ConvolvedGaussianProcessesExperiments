#-------------------------------------------------------------------------------
# SPECIFY DATA
#-------------------------------------------------------------------------------

SOURCE = "Mrk509_2016"


# ❗ Same as experiment 1 but FS set to 400❗

lambda, tobs, yobs, σobs = readdataset(source = SOURCE)


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

FS = 400

include("../../../../runmecommonpart.jl")
