module TradeOff

# First include scripts to make parameters
include("parameters/paras.jl") # Parameters for complete model

using Distributions
using DifferentialEquations
using JLD
using Random

# Include scripts that contain other functions
include("MyPlots.jl") # This means I must always have a version of MyPlots.jl available
include("parameters/make_paras.jl") # Function to make parameters for full model
include("parameters/gen_pool.jl") # Function to generate species pools
include("simulate/sim.jl") # Include simulation code for full model
include("analysis/merge_data.jl") # Allows reconstruction of full trajectories
include("analysis/chemostat.jl") # Allows chemostat testing

# export global constants
export Rgas, ΔGATP, NA

# Declaration of useful constants
global Rgas = 8.31446261815324 # gas constant in J.K^-1.mol^-1
global ΔGATP = 75000.0 # Need to find a reference for this at some point
global NA = 6.02214076e23 # Avogadro's constant in mol−1

# Defining and exporting functions useful in more than one script
export mvector, Keq, Q

# function to generate a vector of values for the maintenance energy requirements m
function mvector(N::Int64, mm::Float64, sdm::Float64)
    @assert mm - 5 * sdm >= 0.0 "This choice could result in negative energy requirements"
    # Initialise vector of m
    m = zeros(N)
    # Make required Gaussian distribution using the provided mean (mm) and SD (sdm)
    d = Normal(mm, sdm)
    for i in 1:N
        m[i] = rand(d)
    end
    return (m)
end

# function to find the equilibrium constant
function Keq(T::Float64, η::Float64, ΔG0::Float64)
    Keq = exp((-ΔG0 - η * ΔGATP) / (Rgas * T))
    return (Keq)
end

# function to find the reaction quotient Q, in the case of 1 to 1 stoichiometry
function Q(S::Float64, P::Float64)
    Q = P / S
    return (Q)
end

# Function to interpolate over a time series
function interpolate_time(ts::Array{Float64,1}, Tg::Float64, T1x::Float64, T2x::Float64)
    return (ts[Tind] * (T1x) / Tg + ts[Tind-1] * (T2x) / Tg)
end

# Function to interpolate over a time series (vectorised form)
function interpolate_time(ts::Array{Float64,2}, Tg::Float64, T1x::Float64, T2x::Float64)
    return (ts[Tind, :] * (T1x) / Tg .+ ts[Tind-1, :] * (T2x) / Tg)
end

end # module
