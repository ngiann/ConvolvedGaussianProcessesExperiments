#-------------------------------------------------------------------------------
# SPECIFY DATA
#-------------------------------------------------------------------------------

# Same as Experiment1 but much more noise and few more observations

SOURCE = "SYNTH"

truemass      = 1e8
trueedfrac    = 20
trueaccretion = TransferFunctions.massaccretionfunction(bhm = truemass, edfrac=trueedfrac, eta=0.1)

lambda, tobs, yobs, σobs, = simulatedatafromgp(mass=truemass, eddingtonfraction=trueedfrac, Tmax = 50, N=50, σ=2.5)

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
TF = pmap(((m,ef),) -> PhysicalTransferFunctionsEddington(mass = m, eddingtonfraction = ef, wavelengths = lambda), Iterators.product(masses, efractions))

FS = 150

include("../../runmecommonpart.jl")
