#-------------------------------------------------------------------------------
# SPECIFY DATA
#-------------------------------------------------------------------------------

# Same as Experiment1 but with higher ρmin = 2.0 and FS = 200
# Also: leave out time series at wavelength A5100 to reduce data size

SOURCE = "Mgc0811"


lambda, tobs, yobs, σobs = readdataset(source = SOURCE)

lambda, tobs, yobs, σobs = lambda[1:end-1], tobs[1:end-1], yobs[1:end-1], σobs[1:end-1]


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

FS = 200

include("../../../runmecommonpart.jl")
