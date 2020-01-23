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

# Function to find K2 from three other rates plus the equilbrium constant, plus checks if reasonable
function findK2(k1::Float64,k2::Float64,K1::Float64,KeQ::Float64,maxr::Float64)
    K2 = (k1*k2)/(K1*KeQ)
    # Test if K2 is a reasonable rate and if not return Inf
    if K2 >= maxr
        return(Inf)
    else
        return(K2)
    end
    return(K2)
end

# function to create alternate parameters from an initial parameter set stored as first row
function shiftks(kset::Array{Float64,2},N::Int64,h::Float64,maxr::Float64)
    # Check if array is appropriate size
    if size(kset,1) != N
        println("Error: Parameter array is too small. Expected $(N) rows and got $(size(kset,1)).")
        error()
    end
    # Generates parameters set at various angles
    for i = 2:N
        # One parameter
        c1 = cos(2*(i-1)*pi/(N-1))
        mc1 = abs(c1)
        dc1 = sign(c1)
        if dc1 == 1.0
            # The 1+mc1*(h-1) keeps the scaling sensible
            kset[i,1] = kset[1,1]*(1+mc1*(h-1))
        else
            kset[i,1] = kset[1,1]/(1+mc1*(h-1))
        end
        # And then the other
        c2 = sin(2*(i-1)*pi/(N-1))
        mc2 = abs(c2)
        dc2 = sign(c2)
        if dc2 == 1.0
            kset[i,3] = kset[1,3]*(1+mc2*(h-1))
        else
            kset[i,3] = kset[1,3]/(1+mc2*(h-1))
        end
        # overwrite rates if greater than maxrate
        if kset[i,1] > maxr
            kset[i,1] = maxr
        end
        if kset[i,3] > maxr
            kset[i,3] = maxr
        end
    end
    return(kset)
end

# function to find the value of KS that maximises
function maxK(k2::Float64,maxr::Float64,KeQ::Float64)
    # Starting values of rates
    k1 = 1.00e3
    K1 = 1.00e3
    # factor to increase/decrease by
    h = 1.05#25
    # Make array to store parameters
    N = 9 # testing 8+1 parameter sets
    kset = fill(Inf,(N,4))
    kset[:,2] .= k2
    kset[1,1] = k1
    kset[1,3] = K1
    KS = fill(Inf,N)
    # Now start while loop
    fsh = false
    while fsh == false
        # Use first row of parameter set to generate other rows
        kset = shiftks(kset,N,h,maxr)
        # Find and add K2 values
        for i = 1:N
            kset[i,4] = findK2(kset[i,1],kset[i,2],kset[i,3],KeQ,maxr)
        end
        # Calculate KS values for the sets
        for i = 1:N
            if isfinite(kset[i,4])
                # Then calculate values for KS
                KS[i] = satK(kset[i,1],kset[i,2],kset[i,3])
            else
                # If K2 value is set Inf then ignore set
                KS[i] = Inf
            end
        end
        # Choose smallest value and make this the new parameter set
        I = argmin(KS)
        # Now overwrite
        kset[1,:] .= kset[I,:]
        # Stop while loop if best parameter set is self
        if I == 1
            fsh = true
        end
    end
    k1 = kset[1,1]
    K1 = kset[1,3]
    K2 = kset[1,4]
    return(k1,K1,K2)
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
    maxr = 1.0e10 # Maximum possible rate
    # Now want to choose maximum populations for a range of k values
    k2s = [1.0,10.0,100.0,1000.0]
    for i = 1:length(k2s)
        k2 = k2s[i]
        k1, K1, K2 = maxK(k2,maxr,KeQ)
        qm, KS, _, kr = qKmY(k1,K1,k2,K2,E0)
        println("Set $(i):")
        println(k2s[i])
        println(η*qm)
        println(m)
        println("S*=$(m*KS/(η*qm-m))")
        QT = Qineq(η,qm,m,kr,KeQ)
        _, _, Ns = stead(KS,kr,η,qm,m,CO,α,δ,θ)
    end
    println(KeQ*0.1)
    return(nothing)
end

@time tradeinvest()
