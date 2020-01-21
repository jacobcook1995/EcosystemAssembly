using Syntrophy
using Plots
using DifferentialEquations
using LaTeXStrings
import PyPlot

# This is a script to write my testing code into
# Anything reusable should be moved into Syntrophy.jl as a seperate callable function

# function to find qm, KS, m, Y based on reference values, kinetic parameters and enzyme concentrations
function qKmY(k1::Float64,K1::Float64,k2::Float64,K2::Float64,E0::Float64)
    # Standard formula for qm, KS, KP and kr
    qm = k2*E0
    KS = (K1+k2)/(k1)
    KP = (K1+k2)/(K2)
    kr = k2/K1
    return(qm,KS,KP,kr)
end

# Function to show threshold for reaction quotient
function Qineq(η::Float64,qm::Float64,m::Float64,kr::Float64,Keq::Float64)
    x = 0.1 # Starting with 10% as a rough guess
    Qi = x*Keq*(η*qm - m)/(η*qm + kr*m)
    return(Qi)
end

# function to find steady state of case without thermodynamic limitation
function stead(KS::Float64,kr::Float64,η::Float64,qm::Float64,m::Float64,CO::Float64,α::Float64,δ::Float64,θ::Float64)
    # Calulate R
    R = m*KS/(η*qm*(1-θ)-m*(1+kr*θ))
    # Then find fractional contribution from S
    S = R/((CO)^6)
    # Then calculate X
    X = (η/m)*(α-δ*S)
    # Then calulate P
    P = 6*m*X/(δ*η)
    return(S,P,X)
end

# Function to do some vague investigating into the tradeoffs
function tradeinvest()
    # First need to define some basic parameters
    # Bunch of assumptions about the environment
    α = 5.55e-6 # These two rates can be changed
    δ = 2.00e-4
    θ = 0.0 # Thermodynamic inhibition not reached
    CO = 0.21 # High value so oxegen isn't limiting
    # Will consider glucose to begin with and then move onto lower free energy changes
    ΔG0 = -2843800.0
    # Assume 38 ATP generated per glucose metabolised
    η = 38.0
    m = 2.16e-19 # maintainance, seems fine shouldn't alter results too much
    Y = 2.36e13 # yield in cells per mole of ATP, scaling constant might need to be redone, but shouldn't alter results
    E0 = 2.5e-20 # Moles of enzyme per cell, not 100% sure of this number, should be improved later
    # Define other thermodynamically relevant constants
    ΔGATP = 75000.0 # Gibbs free energy of formation of ATP in a standard cell
    Temp = 312.0 # Temperature that growth is occuring at in Kelvin
    # Use to find equilbrium constant
    KeQ = Keq(ΔG0,η,ΔGATP,Temp)
    # Now choose k2
    k2 = 140.0 # This is what we had before
    maxr = 1.0e10 # Maximum possible rate
    k1 = maxr
    K1 = 1.00e-5
    K2 = (k1*k2)/(K1*KeQ)
    qm, KS, _, kr = qKmY(k1,K1,k2,K2,E0)
    # Next find threshold
    QT = Qineq(η,qm,m,kr,KeQ)
    # And maximum population
    _, _, Ns = stead(KS,kr,η,qm,m,CO,α,δ,θ)
    println(Ns)
    println(η*α/m)
    return(nothing)
end

@time tradeinvest()
