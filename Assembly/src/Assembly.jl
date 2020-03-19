module Assembly

using Distributions
using DifferentialEquations

# Include scripts that contain other functions
include("MyPlots.jl") # This means I must always have a version of MyPlots.jl available
include("parameters/marsParameters.jl") # Parameters for base case from Marsland et al.
include("parameters/inhibParameters.jl") # Parameters for Inhibition case
include("parameters/make_inhib.jl") # Function to make parameters for model with inhibition
include("simulate/mars.jl") # Include marsland model simulation code
include("simulate/inhib.jl") # Include simulation code for our model

# export global constants
export Rgas, ΔGATP

# Decleration of useful constants
global Rgas = 8.31446261815324 # gas constant in J.K^-1.mol^-1
global ΔGATP = 75000.0 # Need to find a reference for this at some point

# Defining and exporting functions useful in more than one script
export mvector, Keq

# function to generate a vector of values for the maintenance energy requirments m
function mvector(N::Int64,mm::Float64,sdm::Float64)
    @assert mm - 5*sdm >= 0.0 "This choice could result in negative energy requirements"
    # Initialise vector of m
    m = zeros(N)
    # Make required Gaussian distribution using the provided mean (mm) and SD (sdm)
    d = Normal(mm,sdm)
    for i = 1:N
        m[i] = rand(d)
    end
    return(m)
end

# function to find the equilbrium constant
function Keq(T::Float64,η::Float64,ΔG0::Float64)
    Keq = exp((-ΔG0-η*ΔGATP)/(Rgas*T))
    return(Keq)
end

end # module
