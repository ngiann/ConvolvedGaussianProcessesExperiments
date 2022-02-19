#-------------------------------------------------------------------------------
# SPECIFY DATA
#-------------------------------------------------------------------------------

SOURCE = "Ngc2617"


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
TF = [PhysicalTransferFunctionsEddington(mass=m, eddingtonfraction=ef, wavelengths=lambda) for m in masses, ef in efractions]

include("../../../runmecommonpart.jl")
