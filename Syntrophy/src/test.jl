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
    # Need to generate parameter sets with two other rates fixed
    N = 36 # Number of parameter sets
    k1 = zeros(N) # These two can be varied
    K1 = zeros(N)
    K2 = zeros(N) # This one just needed for therodynamic consistency
    # Make parameter sets
    for i = 1:N
        n1 = ceil(i/6)
        k1[i] = 10.0^(2+n1)
        n2 = (i-1) % 6
        K1[i] = 10.0^(n2)
        K2[i] = (k1[i]*k2)/(K1[i]*KeQ)
    end
    # Now make vectors of kinetic parameters
    KS = zeros(N)
    qm = zeros(N)
    kr = zeros(N)
    QT = zeros(N)
    Ns = zeros(N)
    for i = 1:N
        qm[i], KS[i], _, kr[i] = qKmY(k1[i],K1[i],k2,K2[i],E0)
        # Next find threshold
        QT[i] = Qineq(η,qm[i],m,kr[i],KeQ)
        # And maximum population
        _, _, Ns[i] = stead(KS[i],kr[i],η,qm[i],m,CO,α,δ,θ)
    end
    # Remove zero elements
    QT[Ns.<0.0] .= NaN
    Ns[Ns.<0.0] .= NaN
    # Start plotting
    pyplot(dpi=200)
    scatter(QT,Ns)
    savefig("Output/test.png")
    return(nothing)
end

@time tradeinvest()
