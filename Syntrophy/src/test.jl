using Syntrophy
using Plots
using DifferentialEquations
using LaTeXStrings
using SymPy
import PyPlot

# This is a script to write my testing code into
# Anything reusable should be moved into Syntrophy.jl as a seperate callable function

# function to find qm, KS, m, Y based on reference values, kinetic parameters and enzyme concentrations
function qKmY(k1::Float64,K1::Float64,k2::Float64,K2::Float64,E0::Float64)
    # Standard formula for qm, KS, KP and kr
    qm = maxq(k2,E0)
    KS = satK(k1,k2,K1)
    KP = (K1+k2)/(K2)
    kr = krev(k2,K1)
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

# function to find the choice of rates that lead to a maximum QT
function maxQT(k2::Float64,maxrate::Float64,K1K2::Float64,E0::Float64,m::Float64,Keq::Float64,η::Float64)
    # Set number of increments to find
    N = 5000
    # Fraction of maximum QT that is acceptable
    f = 0.001
    # First calculate value for qm (max rate)
    qm = maxq(k2,E0)
    # As starting value use maxrate, make it a logaritmic range
    K2s = exp10.(collect(range(log10(maxrate),stop=log10(K1K2/maxrate),length=N)))
    QT = zeros(N)
    # Loop over each possible K2 value
    for i = 1:N
        # Find corresponding value for K1
        K1 = K1K2/K2s[i]
        # Find reverse rate and use to calculate threshold.
        kr = krev(k2,K1)
        QT[i] = Qineq(η,qm,m,kr,Keq)
    end
    # find maximum value of QT
    mQT = maximum(QT)
    # find first QT within f of it
    I = findfirst(QT .>= mQT-f*mQT)
    # Choose K2 value from the vector
    K2 = K2s[I]
    return(K2)
end

# function to find krates based on max rates, and combined value of K1K2
function krates(maxrate::Float64,K1K2::Float64,k2::Float64,KeQ::Float64,E0::Float64,m::Float64,η::Float64)
    # Find relative value of forward rates
    k1k2 = K1K2*KeQ
    k1 = k1k2/k2
    # Find value of K2 that gives maximum threshold
    K2 = maxQT(k2,maxrate,K1K2,E0,m,KeQ,η)
    K1 = K1K2/K2
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
    # Make vector of Temps to look at
    Temps = collect(Temp-50.0:10.0:Temp+50.0)
    # Setup max rates etc
    maxrate = 1.0e6
    K1K2 = 1.0e7
    # Make vector of k2's to test
    k2 = [collect(0.3:0.1:0.9);collect(1.0:1.0:9.0);collect(10.0:10.0:90.0);collect(100.0:100.0:1000.0)]
    KS = zeros(length(k2),length(Temps))
    qm = zeros(length(k2),length(Temps))
    kr = zeros(length(k2),length(Temps))
    QT = zeros(length(k2),length(Temps))
    Ns = zeros(length(k2),length(Temps))
    KeQ = zeros(length(Temps))
    for j = 1:length(Temps)
        # Calculate equilbrium constant
        KeQ[j] = Keq(ΔG0,η,ΔGATP,Temps[j])
        for i = 1:length(k2)
            # Find k rates from function
            k1, K1, K2 = krates(maxrate,K1K2,k2[i],KeQ[j],E0,m,η)
            qm[i,j], KS[i,j], _, kr[i,j] = qKmY(k1,K1,k2[i],K2,E0)
            QT[i,j] = Qineq(η,qm[i,j],m,kr[i,j],KeQ[j])
            _, _, Ns[i,j] = stead(KS[i,j],kr[i,j],η,qm[i,j],m,CO,α,δ,θ)
        end
    end
    pyplot(dpi=200)
    Lqm = L"q_m"
    plot(qm*1e17,KS,xlabel="$(Lqm) (10^-17)")
    savefig("Output/qmvsKS.png")
    plot(qm*1e17,kr,xlabel="$(Lqm) (10^-17)")
    savefig("Output/qmvskr.png")
    plot(qm*1e17,Ns,xlabel="$(Lqm) (10^-17)")
    savefig("Output/qmvsNs.png")
    plot(qm*1e17,QT,xlabel="$(Lqm) (10^-17)")
    savefig("Output/qmvsQT.png")
    return(nothing)
end

# function to test rate yield tradeoff that Robert suggested
function tradeoff()
    # Define relevant symbols
    N, kp, Δμ0, ΔμATP = symbols("N kp Dm0 DATP")
    # Enter the initial equation
    Q = N*kp*(1-ℯ^(ΔμATP-(Δμ0/N)))
    # Find the differential of this equation
    dQdN = diff(Q,N)
    # Define variables to sub in
    ΔGATP = 75000.0
    ΔG0 = -2843800.0
    k = 1
    # Then sub them in
    sdQdN = subs(dQdN,(kp,k),(ΔμATP,ΔGATP),(Δμ0,-ΔG0))
    # Now plot this for a range of N values
    Ns = collect(30.0:0.01:37.92)
    dQ = zeros(length(Ns))
    for i = 1:length(Ns)
        dQ[i] = subs(sdQdN,N,Ns[i]) |> float
    end
    pyplot(dpi=200)
    plot(Ns,dQ)
    savefig("Output/test.png")
    # Do the same for Q
    sQ = subs(Q,(kp,k),(ΔμATP,ΔGATP),(Δμ0,-ΔG0))
    Qs = zeros(length(Ns))
    for i = 1:length(Ns)
        Qs[i] = subs(sQ,N,Ns[i]) |> float
    end
    plot(Ns,Qs)
    savefig("Output/test2.png")
    # # Then solve for N
    # Nm = SymPy.solve(dQdN,N)
    # println(Nm)
    return(nothing)
end

@time tradeoff()
