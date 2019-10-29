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

# function to return k parameters based on a single k value
function parak(k2::Float64,ΔG0::Float64,η::Float64,ΔGATP::Float64,Temp::Float64)
    K1 = 140.0 # K1 must match value of k2 used for k_R = 1
    # which of k_+1 should I change? k_+1 is the intresting tradeoff
    K2 = 1.28*10.0^(8) # Fix value of K_2
    # Now work out equlibrium constant K in order to find final rate
    K = Keq(ΔG0,η,ΔGATP,Temp)
    k1 = K*K1*K2/(k2)
    return(k1,k2,K1,K2)
end

function maxrate(k2::Float64,ΔG0::Float64,η::Float64,ΔGATP::Float64,Temp::Float64,E0::Float64,Y::Float64,f::Function,u0::Array{Float64,1})
    # Calculate other rates to match k
    k1, k2, K1, K2 = parak(k2,ΔG0,η,ΔGATP,Temp)
    # Use rates to obtain required parameters
    qm, KS, KP, kr = qKmY(k1,K1,k2,K2,E0)
    # Now need to calculate maximal rate mr
    p = [Y,KS,qm,ΔGATP,Temp,kr] # collect parameters
    # put parameters into function
    tspan = (0.0,5000000.0)
    prob = ODEProblem(f,u0,tspan,p)
    # then solve
    sol = solve(prob,adaptive=false,dt=500)
    # Need a check that population has ceased growing
    # Check if pop has changed more than 0.1% in last 10 time steps
    diff = (sol'[end,5]-sol'[end-10,5])/sol'[end,5]
    if diff > 0.001
        println("Not a stable population!")
    end
    # Find maximum population
    mp = maximum(sol'[:,5])
    # Then need to find actual q rate along trajectory
    qs = zeros(length(sol.t))
    for i = 1:length(qs)
        qs[i] = qrate(sol'[i,1:4],KS,qm,ΔGATP,ΔG0,Temp,[-1,-6,6,6],η,kr)
    end
    # Take max of this new vector
    mr = maximum(qs)
    # return both qm and actual maximal observed rate
    return(qm,mr,mp)
end

# function to test whether q_m is a good proxy for maximal growth rate
function testq()
    # Nutrient variables
    α = 5.55*10^(-6)
    δ = 2.00*10^(-4)
    Temp = 312.0 # Temperature that growth is occuring at in Kelvin
    # make nutrients
    # 1 = glucose, 2 = oxegen, 3 = bicarbonate, 4 = hydrogen ion
    nuts = [Nut(1,false,α,δ),Nut(2,true,0,0),Nut(3,false,0,δ),Nut(4,true,0,0)]
    # Now make reactions
    ΔG0 = -2843800.0
    reac = [React(1,[1,2,3,4],[-1,-6,6,6],ΔG0)]
    # physiological value of η
    η = 38.0
    ΔGATP = 75000.0 # Gibbs free energy of formation of ATP in a standard cell
    # Following parameters should not be expected to change between microbes
    E0 = 2.5*10.0^(-20) # Somewhat fudged should be right order of magnitude
    m = 2.16*10.0^(-19) # maintainance
    Y = 2.36*10.0^(13) # yield in cells per mole of ATP
    # Considering 1 microbe with maintaince but no dilution
    mics = Microbe(η,m,1,0.0)
    # Make reduced version of function inputting unchanging microbes
    f(du,u,p,t) = singlepop(du,u,p,nuts,reac,mics,t)
    # Now set initial conditions
    u0 = zeros(5)
    u0[1] = 0.0555 # high initial concentration to ensure growth
    u0[2] = 0.21 # High value so oxegen isn't limiting
    u0[3] = 0.0 # No initial concentration
    u0[4] = 1.00*10.0^(-7) # pH 7
    u0[5] = 100.0
    # Preallocate vectors
    qm = zeros(20)
    mr = zeros(length(qm))
    mp = zeros(length(qm))
    # Loop over vectors
    for i = 1:length(qm)
        # We now want to increase k_{+2} to increase q_m
        k2 = 140.0*i
        qm[i], mr[i], mp[i] = maxrate(k2,ΔG0,η,ΔGATP,Temp,E0,Y,f,u0)
    end
    # Switch backends
    pyplot(dpi=200)
    Ns = L"N^{\ast}"
    p14 = L"10^{14}"
    plot(qm*10.0^17,mp*10.0^(-14),xlabel=L"q_m\;(s^{-1}\,10^{-17})",ylabel="$(Ns) (cells $(p14))")
    savefig("Output/qmvsmp.png")
    plot(qm*10.0^17,mr*10.0^19,xlabel=L"q_m\;(s^{-1}\,10^{-17})",ylabel=L"q\;(s^{-1}\,10^{-19})")
    savefig("Output/qmvsmr.png")
    plot(mr*10.0^19,mp*10.0^(-14),xlabel=L"q\;(s^{-1}\,10^{-19})",ylabel="$(Ns) (cells $(p14))")
    savefig("Output/mrvsmp.png")
    return(nothing)
end

@time testq()
