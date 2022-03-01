#-------------------------------------------------------------------------------
# SPECIFY DATA
#-------------------------------------------------------------------------------

SOURCE = "Mrk509_2017"


lambda, tobs, yobs, σobs = readdataset(source = SOURCE)


#-------------------------------------------------------------------------------
# SPECIFY TRANSFER FUNCTIONS
#-------------------------------------------------------------------------------

# specify kernel
kernelname = "rbf"

# specify physical parameters
masses     = collect(logrange(1e5, 1e10, 64))
efractions = [10.0]

# create combinations of transfer functions
TF = pmap(((m,ef),) -> PhysicalTransferFunctionsEddington(mass = m, eddingtonfraction = ef, wavelengths = lambda), Iterators.product(masses, efractions))

FS = 150
ρmin = 0.1

include("../../../../runmecommonpart.jl")
