
# Purpose of this script is check out a particular model-misspecification scenario:
# This scenario is the case where one of the lightcurves does not bear any connection
# to the others. However, the irrelevant lightcurve has less observations than the others.
# This resembles the MRK509_2016 case.
# We expect the model in this case to set the scaling coefficient of
# the irrelevant lightcurve to zero and still recover the true physical parameters.


#-------------------------------------------------------------------------------
# SPECIFY DATA
#-------------------------------------------------------------------------------

SOURCE = "SYNTH"

truemass      = 1e8
trueedfrac    = 20
trueaccretion = TransferFunctions.massaccretionfunction(bhm = truemass, edfrac=trueedfrac, eta=0.1)

lambda, tobs, yobs, σobs, = simulatedatafromgp(mass=truemass, accretion=trueaccretion, Tmax = 100, N=60, seed = 1)


# generate one more dataset with different seed

lambdaX, tobsX, yobsX, σobsX, = simulatedatafromgp(mass=truemass, accretion=trueaccretion, Tmax = 30, N=20, seed = 10101)

# copy the 4th lightcurve to the original dataset
tobs[4], yobs[4], σobs[4] = tobsX[4], yobsX[4], σobsX[4]

#-------------------------------------------------------------------------------
# SPECIFY TRANSFER FUNCTIONS
#-------------------------------------------------------------------------------

# specify kernel
kernelname = "matern32"

# specify physical parameters
masses     = collect(logrange(1e6, 1e10, 32))
efractions = [10.0; 20.0; 30.0]

# create combinations of transfer functions
TF = pmap(((m,ef),) -> PhysicalTransferFunctionsEddington(mass = m, eddingtonfraction = ef, wavelengths = lambda), Iterators.product(masses, efractions))

include("../../runmecommonpart.jl")
