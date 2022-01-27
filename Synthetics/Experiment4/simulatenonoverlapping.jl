#############################################################################
function simulatenonoverlapping(; seed = 1, mass = 1e6, eddingtonfraction = 10.0, ρ=2.0, σ = 1.0, N = 40,  kernelname="matern32")
#############################################################################

    # σ above scales the σobs defined below

    rg = MersenneTwister(seed)

    #---------------------------------------------------------------------
    # Define array of transfer functions, i.e. one per "band"
    #---------------------------------------------------------------------

    lambda = [4300.0; 5700; 6200]

    tfarray = PhysicalTransferFunctions(mass=mass, eddingtonfraction=eddingtonfraction, wavelengths=lambda)

    scale = [1; 2; 3]
    shift = [5; 6; 7]

    #---------------------------------------------------------------------
    # Specify number of inputs per "band" and sample inputs from two
    # uniform distributions that introduce artificially a gap
    #---------------------------------------------------------------------

    Nobs = [N; N; N] # here we could specify different number of observations for each wavelength

    σobs = [σ*0.10*ones(Nobs[1]),
            σ*0.13*ones(Nobs[2]),
            σ*0.16*ones(Nobs[3])]


    tobs = [rand(Uniform(0.0, 50.0), N), rand(Uniform(50.0, 100.0), N), rand(Uniform(100, 150.0), N)]

    #---------------------------------------------------------------------
    # Define Gaussian process to draw noisy targets
    #---------------------------------------------------------------------

    C = ConvolvedKernelFunction(kernelname, tf = tfarray, ρmax = ρ, fs = 1000)

    setkernelfunction!(C; ρ = ρ, σ² = 1.0)

    K = convolvedScaledCovariance(C, scale,  tobs)

    let

            U, S, V = svd(K)

            K = U * Diagonal(max.(1e-9, abs.(S))) * U'

            makematrixsymmetric!(K)

    end


    #---------------------------------------------------------------------
    # Draw targets and arrange in array
    #---------------------------------------------------------------------

    Y = rand(rg, MvNormal(zeros(sum(Nobs)), K))

    yobs = []

    yobs = [Y[1:Nobs[1]]                                 .+ shift[1] .+ σobs[1].*randn(rg, Nobs[1]),
            Y[1+Nobs[1]:Nobs[1]+Nobs[2]]                 .+ shift[2] .+ σobs[2].*randn(rg, Nobs[2]),
            Y[1+Nobs[1]+Nobs[2]:Nobs[1]+Nobs[2]+Nobs[3]] .+ shift[3] .+ σobs[3].*randn(rg, Nobs[3])]

    figure(-2) ; cla()

    for i in 1:3

      plot(tobs[i], yobs[i], "o", label = @sprintf("%d", lambda[i]))

    end

    legend()

    return tobs, yobs, σobs, tfarray, lambda

end
