#-------------------------------------------------------------------------------
# SPECIFY DATA
#-------------------------------------------------------------------------------

# Purpose of this script is just to verify that we can indeed recover the
# physical parameters (mass and edfraction) in the case of synthetic data.
# Here we control the ground truth and know what the result should be.

SOURCE = "SYNTH"

truemass      = 1e8
trueedfrac    = 20
trueaccretion = TransferFunctions.massaccretionfunction(bhm = truemass, edfrac=trueedfrac, eta=0.1)

lambda, tobs, yobs, σobs, = simulatedatafromgp(mass=truemass, accretion=trueaccretion, Tmax = 100, N=60)

# use only the first three wavelengths
lambda, tobs, yobs, σobs = lambda[1:3], tobs[1:3], yobs[1:3], σobs[1:3]

#-------------------------------------------------------------------------------
# SPECIFY TRANSFER FUNCTIONS
#-------------------------------------------------------------------------------

# specify kernel
kernelname = "matern32"

# specify physical parameters
masses     = collect(logrange(1e6, 1e10, 32))
efractions = [10.0; 20.0; 30.0]


# create combinations of transfer functions
TF = pmap(((m,ef),) -> ShakuraBoxTransferFunctionsEddington(mass = m, eddingtonfraction = ef, wavelengths = lambda), Iterators.product(masses, efractions))

FS = 200

ρmin = 1.0

include("../../runmecommonpart.jl")
