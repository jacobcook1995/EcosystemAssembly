using Syntrophy
using Plots
using DifferentialEquations
using LaTeXStrings
# Not needed at moment
# import PyPlot

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
    # p[2] = KS, p[3] = qm, p[4] = ΔGATP, p[5] = Temp
    q = qrate(u[1:N],p[2],p[3],p[4],ΔG0,p[5],stc,η)
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

function gluc()
    # Nutrient variables
    α = 5.55*10^(-6)
    δ = 2.00*10^(-4) # 1.00*10^(-4)
    # make nutrients
    # 1 = glucose, 2 = oxegen, 3 = bicarbonate, 4 = hydrogen ion
    nuts = [Nut(1,false,α,δ),Nut(2,true,0,0),Nut(3,false,0,δ),Nut(4,true,0,0)]
    # Now make reactions
    ΔG0 = -2843800.0
    reac = [React(1,[1,2,3,4],[-1,-6,6,6],ΔG0)]
    # microbe variables
    η = 15.0
    r = 1 # Only reaction
    m = 2.16*10^(-19) # maintainance
    # Considering 1 microbe with no maintaince and no dilution
    mics = Microbe(η,m,1,0.0)
    # Set intial populations and nutrient concentrations
    pops = 100.0
    concs = zeros(length(nuts))
    # define initial concentrations
    concs[1] = 0.0555 # high initial concentration to ensure growth
    concs[2] = 0.21 # High value so oxegen isn't limiting
    concs[3] = 0.0 # No initial concentration
    concs[4] = 1.00*10.0^(-7) # pH 7
    # Define some constants
    ΔGATP = 75000.0 # Gibbs free energy of formation of ATP in a standard cell
    Temp = 312.0 # Temperature that growth is occuring at in Kelvin
    # All other constants require quite a bit of defining
    Y = 2.36*10.0^(13) # yield in cells per mole of ATP
    KS = 2.40*10.0^(-5) # saturation constant (substrate)
    qm = 3.42*10.0^(-18) # maximal rate substrate consumption mol cell s^-1
    p = [Y,KS,qm,ΔGATP,Temp]
    u0 = [concs;pops]
    tspan = (0.0,5000000.0)

    # Make reduced version of function inputting unchanging microbes
    f(du,u,p,t) = singlepop(du,u,p,nuts,reac,mics,t)
    prob = ODEProblem(f,u0,tspan,p)
    sol = solve(prob,adaptive=false,dt=500) # turned dt down to make plots look nicer
    # Now do plotting
    plot(sol.t,sol'[:,1])
    savefig("Output/$(η)test1.png")
    plot(sol.t,sol'[:,3])
    savefig("Output/$(η)test3.png")
    plot(sol.t,sol'[:,5])
    savefig("Output/$(η)test5.png")
    stoc = (reac.↦:stc)[1]
    # Print final thermodynamic term
    L = size(sol',1)
    θs = zeros(L)
    for i = 1:L
        θs[i] = θT(sol'[i,1:4],stoc,ΔGATP,ΔG0,η,Temp)
    end
    plot(sol.t,θs)
    savefig("Output/Theta.png")
    return(nothing)
end

# function to find qm, KS, m, Y based on reference values, kinetic parameters and enzyme concentrations
function qKmY(k1::Float64,K1::Float64,k2::Float64,K2::Float64,E0::Float64,E0ref::Float64,mref::Float64,Yref::Float64)
    # Standard formula for qm and KS
    qm = k2*E0
    KS = (K1+k2)/(k1)
    # poportional sub-linear increase in maintainance with amount
    m = mref + 0.1*mref*(E0-E0ref)/E0ref # More than a 10% penalty and the tradeoff is lost
    # decrease in same proportion for yield
    Y = Yref - 0.1*Yref*(E0-E0ref)/E0ref
    return(qm,KS,m,Y)
end

# function to investigate the r vs K tradeoff to attempt to see if it has a thermodynamic origin
function rvsK()
    # Set up for glucose respiration => A lot of these things need changing
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
    # microbe variables
    η = 37.5
    ΔGATP = 75000.0 # Gibbs free energy of formation of ATP in a standard cell
    r = 1 # Only reaction
    # Define kinetic parameters explicitly
    E0ref = 2.5*10.0^(-20) # Somewhat fudged should be right order of magnitude
    k2 = 140.0 # Infered from old qm
    K1 = k2 # As we assume that both directions have same maximal rate
    k1 = 1.17*10.0^(7) # Infered from KS based on value of the above
    # Now work out equlibrium constant K in order to find final rate
    K = Keq(ΔG0,η,ΔGATP,Temp)
    K2 = k1*k2/(K1*K)
    # These can be used as reference values corresponding to the case of E0ref
    mr = 2.16*10.0^(-19) # maintainance
    Yr = 2.36*10.0^(13) # yield in cells per mole of ATP
    # Set intial populations and nutrient concentrations
    pops = 100.0
    concs = zeros(length(nuts))
    # define initial concentrations
    concs[1] = 0.0555 # high initial concentration to ensure growth
    concs[2] = 0.21 # High value so oxegen isn't limiting
    concs[3] = 0.0 # No initial concentration
    concs[4] = 1.00*10.0^(-7) # pH 7
    u0 = [concs;pops]
    tspan = (0.0,5000000.0)
    # Important difference now is in the value of E0
    E0 = E0ref
    # Find KS, qm, maintainance and yield using function
    qm, KS, m, Y = qKmY(k1,K1,k2,K2,E0,E0ref,mr,Yr)
    # Considering 1 microbe with maintaince but no dilution
    mics = Microbe(η,m,1,0.0)
    p = [Y,KS,qm,ΔGATP,Temp]
    # Make reduced version of function inputting unchanging microbes
    f(du,u,p,t) = singlepop(du,u,p,nuts,reac,mics,t)
    prob = ODEProblem(f,u0,tspan,p)
    sol = solve(prob,adaptive=false,dt=500) # turned dt down to make plots look nicer
    println("K = $(maximum(sol'[:,5]))")
    # Now find maximal growth rate r
    v = zeros(length(sol'[:,5])-1)
    for i = 2:length(sol'[:,5])
        v[i-1] = (sol'[i,5]-sol'[i-1,5])/500
    end
    println("r = $(maximum(v))")
    # Important difference now is in the value of E0
    E0 = 2*E0ref
    # Find KS, qm, maintainance and yield using function
    qm, KS, m, Y = qKmY(k1,K1,k2,K2,E0,E0ref,mr,Yr)
    # Considering 1 microbe with maintaince but no dilution
    mics = Microbe(η,m,1,0.0)
    p = [Y,KS,qm,ΔGATP,Temp]
    # Make reduced version of function inputting unchanging microbes
    f(du,u,p,t) = singlepop(du,u,p,nuts,reac,mics,t)
    prob = ODEProblem(f,u0,tspan,p)
    sol = solve(prob,adaptive=false,dt=500) # turned dt down to make plots look nicer
    println("K = $(maximum(sol'[:,5]))")
    for i = 2:length(sol'[:,5])
        v[i-1] = (sol'[i,5]-sol'[i-1,5])/500
    end
    println("r = $(maximum(v))")
    return(nothing)
end

# function to return k parameters based on a single k value
function parak(k1::Float64,k2::Float64,ΔG0::Float64,η::Float64,ΔGATP::Float64,Temp::Float64)
    K1 = k2 # As we assume that both directions have same maximal rate
    # Now work out equlibrium constant K in order to find final rate
    K = Keq(ΔG0,η,ΔGATP,Temp)
    K2 = k1*k2/(K1*K)
    return(k1,k2,K1,K2)
end

# a second similar function to test the effect of changing reaction rates
function rvsK2()
    # Set up for glucose respiration => A lot of these things need changing
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
    # microbe variables
    η = 37.5
    ΔGATP = 75000.0 # Gibbs free energy of formation of ATP in a standard cell
    r = 1 # Only reaction
    # Define kinetic parameters explicitly
    E0 = 2.5*10.0^(-20) # Somewhat fudged should be right order of magnitude
    # These can be used as reference values corresponding to the case of E0ref
    mr = 2.16*10.0^(-19) # maintainance
    Yr = 2.36*10.0^(13) # yield in cells per mole of ATP
    # Set intial populations and nutrient concentrations
    pops = 100.0
    concs = zeros(length(nuts))
    # define initial concentrations
    concs[1] = 0.0555 # high initial concentration to ensure growth
    concs[2] = 0.21 # High value so oxegen isn't limiting
    concs[3] = 0.0 # No initial concentration
    concs[4] = 1.00*10.0^(-7) # pH 7
    u0 = [concs;pops]
    tspan = (0.0,5000000.0)
    # The important difference now is in the value of k1
    k1 = 1.17*10.0^(7) # Set default k1 here
    k2 = 140.0
    k1, k2, K1, K2 = parak(k1,k2,ΔG0,η,ΔGATP,Temp)
    println("k1 = $(k1)")
    println("k2 = $(k2)")
    println("K1 = $(K1)")
    println("K2 = $(K2)")
    # Find KS, qm, maintainance and yield using function
    qm, KS, m, Y = qKmY(k1,K1,k2,K2,E0,E0,mr,Yr)
    # Considering 1 microbe with maintaince but no dilution
    mics = Microbe(η,m,1,0.0)
    p = [Y,KS,qm,ΔGATP,Temp]
    println("KS = $(KS)")
    # Make reduced version of function inputting unchanging microbes
    f(du,u,p,t) = singlepop(du,u,p,nuts,reac,mics,t)
    prob = ODEProblem(f,u0,tspan,p)
    sol = solve(prob,adaptive=false,dt=500) # turned dt down to make plots look nicer
    println("K = $(maximum(sol'[:,5]))")
    # Now find maximal growth rate r
    v = zeros(length(sol'[:,5])-1)
    for i = 2:length(sol'[:,5])
        v[i-1] = (sol'[i,5]-sol'[i-1,5])/500
    end
    println("r = $(maximum(v))")
    # Important difference now is in the value of E0
    k1 = (1/4)*k1
    k1, k2, K1, K2 = parak(k1,k2,ΔG0,η,ΔGATP,Temp)
    println("k1 = $(k1)")
    println("k2 = $(k2)")
    println("K1 = $(K1)")
    println("K2 = $(K2)")
    # Find KS, qm, maintainance and yield using function
    qm, KS, m, Y = qKmY(k1,K1,k2,K2,E0,E0,mr,Yr)
    # Considering 1 microbe with maintaince but no dilution
    mics = Microbe(η,m,1,0.0)
    p = [Y,KS,qm,ΔGATP,Temp]
    println("KS = $(KS)")
    # Make reduced version of function inputting unchanging microbes
    f(du,u,p,t) = singlepop(du,u,p,nuts,reac,mics,t)
    prob = ODEProblem(f,u0,tspan,p)
    sol = solve(prob,adaptive=false,dt=500) # turned dt down to make plots look nicer
    println("K = $(maximum(sol'[:,5]))")
    for i = 2:length(sol'[:,5])
        v[i-1] = (sol'[i,5]-sol'[i-1,5])/500
    end
    println("r = $(maximum(v))")
    plot(sol'[:,5])
    savefig("Output/test1.png")
    plot(sol'[:,3])
    savefig("Output/test2.png")
    plot(sol'[:,1])
    savefig("Output/test3.png")
    return(nothing)
end

@time rvsK2()
