using Syntrophy
using Plots
using DifferentialEquations
using LaTeXStrings
import PyPlot

# This is a script to write my testing code into
# Anything reusable should be moved into Syntrophy.jl as a seperate callable function

# Function to update population and nutrient concentrations
# This is run for a single population utilising a single reaction
function singlepop(du::Array{Float64,1},u::Array{Float64,1},p::Array{Float64,1},nuts::Array{Nut,1},reacs::Array{React,1},
                mic::Microbe,t::Float64)
    # Extract required parameters
    Y = p[1]
    # Extract relevant data from vector of nutrients
    α = nuts.↦:α
    δ = nuts.↦:δ
    con = nuts.↦:cst
    N = length(nuts) # Number of nutrients
    # And relevant data from vector of microbes
    η = mic.η
    m = mic.m # running for single microbe
    M = 1 # Number of microbes
    # Extract reaction stochiometry
    stc = (reacs.↦:stc)[1]
    ΔG0 = (reacs.↦:ΔG0)[1]
    # Now calculate q
    # p[2] = KS, p[3] = qm, p[4] = ΔGATP, p[5] = Temp, p[6] = kr
    q = qrate(u[1:N],p[2],p[3],p[4],ΔG0,p[5],stc,η,p[6])
    # Make vector to store nutrient changes due to consumption
    δX = zeros(N)
    for i = 1:length(stc)
        for j = N+1:N+M
            δX[i] += stc[i]*q*u[j] # Assumes single reaction
        end
    end
    # q has no dependance on population
    # Now update nutrients
    for i = 1:N
        if con[i] == false
            du[i] = α[i]-δ[i]*u[i]+δX[i]
        else
            du[i] = 0
        end
    end
    # Then calculate population changes
    for i = N+1:N+M
        j = i-N
        E = netE(η,q,m)
        if E >= 0.0 # find if growing or decaying
            du[i] = E*Y*u[i] # No dilution rate so can ignore
        else
            du[i] = E*Y*u[i]
        end
    end
    return(du)
end

# function to find qm, KS, m, Y based on reference values, kinetic parameters and enzyme concentrations
function qKmY(k1::Float64,K1::Float64,k2::Float64,K2::Float64,E0::Float64)
    # Standard formula for qm, KS, KP and kr
    qm = k2*E0
    KS = (K1+k2)/(k1)
    KP = (K1+k2)/(K2)
    kr = k2/K1
    return(qm,KS,KP,kr)
end

# Function to calulate the value that needs to be exceeded for thermodynamic effects
# This is stored here as it is a function that has a mainly illustrative purpose
function Qineq(η::Float64,qm::Float64,m::Float64,kr::Float64,Keq::Float64)
    x = 0.1 # Starting with 10% as a rough guess
    Qi = x*Keq*(η*qm - m)/(η*qm + kr*m)
    return(Qi)
end

# Maybe just use a function to find threshold θ instead if we're just intrested in showing the tradeoff
function θineq(η::Float64,qm::Float64,m::Float64,kr::Float64)
    x = 0.1 # Starting with 10% as a rough guess
    θi = x*(η*qm - m)/(η*qm + kr*m)
    return(θi)
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

# Increase all rates simulatanously
function rateincrease(set::Array{Float64,2})
    # Find length of parameter set
    L = size(set,1)
    # check that the initial parameter set is defined
    if set[1,1] <= 0.0 || set[1,2] <= 0.0 || set[1,3] <= 0.0 || set[1,4] <= 0.0
        println("Error: Initial parameter set provided incorrectly.")
    end
    # Now fill in remainder of parameter set
    for i = 2:L
        set[i,1] = i*set[1,1]
        set[i,2] = i*set[1,2]
        set[i,3] = i*set[1,3]
        set[i,4] = i*set[1,4]
    end
    return(set)
end

# Increase rate of S binding and unbinding
function step1increase(set::Array{Float64,2})
    # Find length of parameter set
    L = size(set,1)
    # check that the initial parameter set is defined
    if set[1,1] <= 0.0 || set[1,2] <= 0.0 || set[1,3] <= 0.0 || set[1,4] <= 0.0
        println("Error: Initial parameter set provided incorrectly.")
    end
    # Now fill in remainder of parameter set
    for i = 2:L
        set[i,1] = i*set[1,1]
        set[i,2] = set[1,2]
        set[i,3] = i*set[1,3]
        set[i,4] = set[1,4]
    end
    return(set)
end

# Increase rate of P unbinding and binding
function step2increase(set::Array{Float64,2})
    # Find length of parameter set
    L = size(set,1)
    # check that the initial parameter set is defined
    if set[1,1] <= 0.0 || set[1,2] <= 0.0 || set[1,3] <= 0.0 || set[1,4] <= 0.0
        println("Error: Initial parameter set provided incorrectly.")
    end
    # Now fill in remainder of parameter set
    for i = 2:L
        set[i,1] = set[1,1]
        set[i,2] = i*set[1,2]
        set[i,3] = set[1,3]
        set[i,4] = i*set[1,4]
    end
    return(set)
end

# Increase both binding rates
function bindincrease(set::Array{Float64,2})
    # Find length of parameter set
    L = size(set,1)
    # check that the initial parameter set is defined
    if set[1,1] <= 0.0 || set[1,2] <= 0.0 || set[1,3] <= 0.0 || set[1,4] <= 0.0
        println("Error: Initial parameter set provided incorrectly.")
    end
    # Now fill in remainder of parameter set
    for i = 2:L
        set[i,1] = i*set[1,1]
        set[i,2] = set[1,2]
        set[i,3] = set[1,3]
        set[i,4] = i*set[1,4]
    end
    return(set)
end

# Increase both unbinding rates
function unbindincrease(set::Array{Float64,2})
    # Find length of parameter set
    L = size(set,1)
    # check that the initial parameter set is defined
    if set[1,1] <= 0.0 || set[1,2] <= 0.0 || set[1,3] <= 0.0 || set[1,4] <= 0.0
        println("Error: Initial parameter set provided incorrectly.")
    end
    # Now fill in remainder of parameter set
    for i = 2:L
        set[i,1] = set[1,1]
        set[i,2] = i*set[1,2]
        set[i,3] = i*set[1,3]
        set[i,4] = set[1,4]
    end
    return(set)
end

# One of two actual tradeoffs, increase substrate binding decrease product unbinding
function rvsKtrade(set::Array{Float64,2})
    # Find length of parameter set
    L = size(set,1)
    # check that the initial parameter set is defined
    if set[1,1] <= 0.0 || set[1,2] <= 0.0 || set[1,3] <= 0.0 || set[1,4] <= 0.0
        println("Error: Initial parameter set provided incorrectly.")
    end
    # Now fill in remainder of parameter set
    for i = 2:L
        set[i,1] = i*set[1,1]
        set[i,2] = set[1,2]/i
        set[i,3] = set[1,3]
        set[i,4] = set[1,4]
    end
    return(set)
end

# Other actual tradeoff, increase substrate unbinding decrease product binding
function Kvsθtrade(set::Array{Float64,2})
    # Find length of parameter set
    L = size(set,1)
    # check that the initial parameter set is defined
    if set[1,1] <= 0.0 || set[1,2] <= 0.0 || set[1,3] <= 0.0 || set[1,4] <= 0.0
        println("Error: Initial parameter set provided incorrectly.")
    end
    # Now fill in remainder of parameter set
    for i = 2:L
        set[i,1] = set[1,1]
        set[i,2] = set[1,2]
        set[i,3] = i*set[1,3]
        set[i,4] = set[1,4]/i
    end
    return(set)
end

# function to find maximum rate r, maximum population N and threshold for thermodynamic inhibition θt
function maximas(set::Array{Float64,2},E0::Float64,η::Float64,m::Float64,CO::Float64,α::Float64,δ::Float64)
    # Find length of parameter set
    L = size(set,1)
    # make vectors to store quantities
    r = zeros(L)
    KS = zeros(L)
    kr = zeros(L)
    N = zeros(L)
    θt = zeros(L)
    # Now find the various two maxima and the trheshold
    for i = 1:L
        # Normal method of finding parameters
        r[i], KS[i], _, kr[i] = qKmY(set[i,1],set[i,3],set[i,2],set[i,4],E0)
        # θ = 0 as want max without thermodynamic inhibition
        _, _, N[i] = stead(KS[i],kr[i],η,r[i],m,CO,α,δ,0.0)
        # Now find threhold for thermodynamic inhibition
        θt[i] = θineq(η,r[i],m,kr[i])
    end
    return(r,N,θt)
end

# Function to demonstrate all 7 tradeoffs
function tradeoffs()
    # This is the parameterisation I used initially want to check that these values of q_m and KS are the same
    # Nutrient variables
    α = 5.55e-6
    δ = 2.00e-4
    # make nutrients
    # 1 = glucose, 2 = oxegen, 3 = bicarbonate, 4 = hydrogen ion
    nuts = [Nut(1,false,α,δ),Nut(2,true,0,0),Nut(3,false,0,δ),Nut(4,true,0,0)]
    # Now make reactions
    ΔG0 = -2843800.0
    reac = [React(1,[1,2,3,4],[-1,-6,6,6],ΔG0)]
    # microbe variables
    η = 38.0
    r = 1 # Only reaction
    m = 2.16e-19 # maintainance
    Y = 2.36e13 # yield in cells per mole of ATP
    E0 = 2.5e-20 # Somewhat fudged should be right order of magnitude
    # Considering 1 microbe with no maintaince and no dilution
    mics = [Microbe(η,m,r,0.0)]
    # Set intial populations and nutrient concentrations
    u0 = zeros(length(nuts)+length(mics)) # make vector of initial conditions
    # define initial concentrations
    u0[1] = 0.0555 # high initial concentration to ensure growth
    u0[2] = 0.21 # High value so oxegen isn't limiting
    u0[3] = 0.0 # No initial concentration
    u0[4] = 1.00e-7 # pH 7
    u0[5] = 100.0 # 100 hundred initial cells
    # Define other thermodynamically relevant constants
    ΔGATP = 75000.0 # Gibbs free energy of formation of ATP in a standard cell
    Temp = 312.0 # Temperature that growth is occuring at in Kelvin
    # New parameters for base case
    K1 = 140.0
    k2 = 140.0
    k1 = 1.17e7
    # Defining three above parameters fix q_m, k_r and K_S
    # Now set the 4th rate K2 (and thus K_P) to ensure thermodynamic consistency
    KeQ = Keq(ΔG0,η,ΔGATP,Temp)
    K2 = k1*k2/(KeQ*K1)
    # Group these initial parameters together
    ki = [k1,k2,K1,K2]
    # Now want to define 7 types of changes
    ps = zeros(7,20,4)
    # Asign these parameters as first parameters for each tradeoff
    for i = 1:7
        ps[i,1,1:4] = ki
    end
    # Now fill out parameter sets by providing to a function for each tradeoff
    ps[1,:,:] = rateincrease(ps[1,:,:])
    ps[2,:,:] = step1increase(ps[2,:,:])
    ps[3,:,:] = step2increase(ps[3,:,:])
    ps[4,:,:] = bindincrease(ps[4,:,:])
    ps[5,:,:] = unbindincrease(ps[5,:,:])
    ps[6,:,:] = rvsKtrade(ps[6,:,:])
    ps[7,:,:] = Kvsθtrade(ps[7,:,:])
    # For now just want data on max growth rate, max population (without thermodynamic inhibition),
    # and threshold for thermodynamic inhibition
    r = zeros(size(ps,1),size(ps,2)) # maximal rate
    N = zeros(size(ps,1),size(ps,2)) # maximal populations
    θt = zeros(size(ps,1),size(ps,2)) # threshold θ for thermodynamic inhibition
    for i = 1:size(ps,1)
        r[i,:], N[i,:], θt[i,:] = maximas(ps[i,:,:],E0,η,m,u0[2],α,δ)
    end
    # Now plot these tradeoffs so that they can be seen
    pyplot(dpi=200)
    plot(r[6,:]*1.0e17,N[6,:])
    savefig("Output/rvsK.png")
    plot(r[6,:]*1.0e17,θt[6,:])
    savefig("Output/rvsKT.png")
    plot(N[7,:],r[7,:]*1.0e18)
    savefig("Output/θvsKrate.png")
    plot(N[7,:],θt[7,:])
    savefig("Output/θvsK.png")
    return(nothing)
end
# Then function that takes parameter sets for each tradeoff and simulates and makes graphs of tradeoff
# Then loop over for the seven tradeoffs

@time tradeoffs()
