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
    KS = satK(k1,k2,K1)
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

# function to find krates based on max rates, and combined value of K1K2
function krates(maxrate::Float64,K1K2::Float64,k2::Float64,KeQ::Float64)
    # Find relative value of forward rates
    k1k2 = K1K2*KeQ
    k1 = k1k2/k2
    # maximise one reverse rate
    K1 = maxrate
    K2 = K1K2/K1
    return(k1,K1,K2)
end

# Function to test simple idea
function fixK1K2()
    # First need to define some basic parameters
    # Bunch of assumptions about the environment
    α = 5.55e-4 # These two rates can be changed
    δ = 1.00e-4
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
    # Setup max rates etc
    maxrate = 1.0e6
    K1K2 = 1.0e8
    # Make vector of k2's to test
    k2 = collect(10.0:10.0:10000.0)
    k1 = zeros(length(k2))
    K1 = zeros(length(k2))
    K2 = zeros(length(k2))
    KS = zeros(length(k2))
    qm = zeros(length(k2))
    kr = zeros(length(k2))
    QT = zeros(length(k2))
    Ns = zeros(length(k2))
    for i = 1:length(k2)
        # Find k rates from function
        k1[i], K1[i], K2[i] = krates(maxrate,K1K2,k2[i],KeQ)
        qm[i], KS[i], _, kr[i] = qKmY(k1[i],K1[i],k2[i],K2[i],E0)
        QT[i] = Qineq(η,qm[i],m,kr[i],KeQ)
        _, _, Ns[i] = stead(KS[i],kr[i],η,qm[i],m,CO,α,δ,θ)
    end
    pyplot(dpi=200)
    Lqm = L"q_m"
    plot(qm*1e18,KS,xlabel="$(Lqm) (10^-18)")
    savefig("Output/qmvsKS.png")
    plot(qm*1e18,kr,xlabel="$(Lqm) (10^-18)")
    savefig("Output/qmvskr.png")
    plot(qm*1e18,Ns,xlabel="$(Lqm) (10^-18)")
    savefig("Output/qmvsNs.png")
    plot(qm*1e18,QT,xlabel="$(Lqm) (10^-18)")
    savefig("Output/qmvsQT.png")
    return(nothing)
end

@time fixK1K2()
