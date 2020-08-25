module Assembly

# First include scripts to make parameters
include("parameters/inhibParameters.jl") # Parameters for Inhibition case
include("parameters/protParameters.jl") # Parameters for proteome case
include("parameters/fullParameters.jl") # Parameters for complete model

using Distributions
using DifferentialEquations
using SymPy
using LinearAlgebra

# Include scripts that contain other functions
include("MyPlots.jl") # This means I must always have a version of MyPlots.jl available
include("parameters/make_inhib.jl") # Function to make parameters for model with inhibition
include("parameters/make_prot.jl") # Function to make parameters for proteome model
include("parameters/make_full.jl") # Function to make parameters for full model
include("simulate/inhib.jl") # Include simulation code for our model
include("simulate/analytic.jl") # Include analytic functions to assit simulations
include("simulate/proteome.jl") # Include simulation code for proteome model
include("simulate/full.jl") # Include simulation code for full model

# export global constants
export Rgas, ΔGATP, NA

# Decleration of useful constants
global Rgas = 8.31446261815324 # gas constant in J.K^-1.mol^-1
global ΔGATP = 75000.0 # Need to find a reference for this at some point
global NA = 6.02214076e23 # Avogadro's constants in mol−1

# Defining and exporting functions useful in more than one script
export mvector, Keq, Q

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

# function to find the reaction quotient Q, in the case of 1 to 1 stochiometery
function Q(S::Float64,P::Float64)
    Q = P/S
    return(Q)
end

end # module
