#-------------------------------------------------------------------------------
# SPECIFY DATA
#-------------------------------------------------------------------------------

SOURCE = "Mrk509_2017"


lambda, tobs, yobs, σobs = readdataset(source = SOURCE)


#-------------------------------------------------------------------------------
# SPECIFY TRANSFER FUNCTIONS
#-------------------------------------------------------------------------------

# specify kernel
kernelname = "matern32"

# specify physical parameters
masses     = collect(logrange(1e5, 6e10, 64))
efractions = [1;5;10.0]

# create combinations of transfer functions
TF = @showprogress [PhysicalTransferFunctionsEddington(mass=m, eddingtonfraction=ef, wavelengths=lambda) for m in masses, ef in efractions]

include("../../../../runmecommonpart.jl")
