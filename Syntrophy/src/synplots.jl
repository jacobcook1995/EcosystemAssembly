using Syntrophy
using Plots
using DifferentialEquations
using LaTeXStrings
import PyPlot

# This is a script to store the functions that plot the figures shown in my Syntrophy SI

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

# A function to make a plot of the maxcosump rate for a microbial population in glucose
function maxcosump()
    # Nutrient variables
    α = 5.55*10^(-6)
    δ = 2.00*10^(-4) #1.00*10^(-4)
    # make nutrients
    # 1 = glucose, 2 = oxegen, 3 = bicarbonate, 4 = hydrogen ion
    nuts = [Nut(1,false,α,δ),Nut(2,true,0,0),Nut(3,false,0,δ),Nut(4,true,0,0)]
    # Now make reactions
    ΔG0 = -2843800.0
    reac = [React(1,[1,2,3,4],[-1,-6,6,6],ΔG0)]
    # make set of microbes
    ηt = [collect(0:0.5:3.5);collect(4:4:32);collect(36:39);collect(39.25:0.25:43)]
    N = length(ηt) # (maximal η)-1
    r = 1 # Only reaction
    m = 2.16*10^(-19) # maintainance
    mics = Array{Microbe,1}(undef,N)
    for i = 1:N
        mics[i] = Microbe(ηt[i],m,r,0.0)
    end
    # Set intial population and nutrient concentrations for all cases
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
    u0 = [concs;pops] # u0 and p same for every microbe
    tspan = (0.0,5000000.0)
    stoc = (reac.↦:stc)[1]

    # Preallocate vector and calculate maxium ATP generation rate
    atpgen = zeros(N)
    ηs = zeros(N)
    θs = zeros(N)
    for i = 1:N
        settle = false
        # New version of function each time for each microbe
        f(du,u,p,t) = singlepop(du,u,p,nuts,reac,mics[i],t)
        # This step will take a while
        prob = ODEProblem(f,u0,tspan,p)
        sol = solve(prob,adaptive=false,dt=100) # turned dt down to make plots look nicer
        # Check if pop has changed more than 0.1% in last 10 time steps
        diff = (sol'[end,5]-sol'[end-10,5])/sol'[end,5]
        if diff < 0.001
            settle = true
        end
        # Print if it hasn't settled down to sufficently steady popultaion
        if settle == false
            println("Problem with eta = $(ηt[i])")
            println(diff)
        end
        # Now calulate max ATP genration rate
        ηc = mics[i].η # ATP usage amount
        qr = qrate(sol'[end,1:4],KS,qm,ΔGATP,ΔG0,Temp,stoc,ηc) # reaction rate
        atpgen[i] = sol'[end,5]*ηc*qr # multiply by population to get total rate
        θs[i] = θT(sol'[end,1:4],stoc,ΔGATP,ΔG0,ηc,Temp) # obtain thermodynamic term
        ηs[i] = ηc # save η for later plotting
    end
    pyplot(dpi=150)
    # Plot maximal ATP generation rate
    plot(ηs,atpgen,xaxis=L"\eta\;\;mol_{ATP}\;(mol_{reaction})^{-1}",yaxis="Maximal ATP production rate mol/s",label="")
    plot!(title="ATP rate vs efficency trade off")
    savefig("Output/Atpgenrate.png")
    # Plot thermodynamicinhibition term
    plot(ηs,θs,xaxis=L"\eta\;\;mol_{ATP}\;(mol_{reaction})^{-1}",yaxis=L"\theta\;\;",label="")
    plot!(title="Thermodynamic inhibition as a function of efficency")
    savefig("Output/ThermInhib.png")
    return(nothing)
end

function limunlim()
    # Nutrient variables
    α = 5.55*10^(-6)
    δ = 2.00*10^(-4) #1.00*10^(-4)
    # make nutrients
    # 1 = glucose, 2 = oxegen, 3 = bicarbonate, 4 = hydrogen ion
    nuts = [Nut(1,false,α,δ),Nut(2,true,0,0),Nut(3,false,0,δ),Nut(4,true,0,0)]
    # Now make reactions
    ΔG0 = -2843800.0
    reac = [React(1,[1,2,3,4],[-1,-6,6,6],ΔG0)]
    # microbe variables
    ηs = [32.0,41.5]
    r = 1 # Only reaction
    m = 2.16*10^(-19) # maintainance
    # Considering 1 microbe with no maintaince and no dilution
    mics = [Microbe(ηs[1],m,r,0.0),Microbe(ηs[2],m,r,0.0)]
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
    stoc = (reac.↦:stc)[1]

    # Make reduced version of function inputting unchanging microbes
    # Do non-limited case first
    pyplot(dpi=150,size =(500,500))
    f(du,u,p,t) = singlepop(du,u,p,nuts,reac,mics[1],t)
    prob = ODEProblem(f,u0,tspan,p)
    sol = solve(prob,adaptive=false,dt=500) # turned dt down to make plots look nicer
    # Plot nutrient concentrations on same graph
    p1 = plot(sol.t,[sol'[:,1],sol'[:,3]],ylabel="Concentration M L^-1",label=["substrate" "product"],title="No thermodynamic inhibition")
    # Then plot population on another subplot
    p2 = plot(sol.t,sol'[:,5],label="",ylabel="population cells L^-1")
    plot(p1,p2,layout=(2,1),xlabel="time s")
    savefig("Output/Unlim.png")

    # Then do same for limited case
    g(du,u,p,t) = singlepop(du,u,p,nuts,reac,mics[2],t)
    prob = ODEProblem(g,u0,tspan,p)
    sol = solve(prob,adaptive=false,dt=500) # turned dt down to make plots look nicer
    # Plot nutrient concentrations on same graph
    p1 = plot(sol.t,[sol'[:,1],sol'[:,3]],ylabel="Concentration M L^-1",label=["substrate" "product"],title="Thermodynamic inhibition")
    # Then plot population on another subplot
    p2 = plot(sol.t,sol'[:,5],label="",ylabel="population cells L^-1")
    plot(p1,p2,layout=(2,1),xlabel="time s")
    savefig("Output/Lim.png")
    return(nothing)
end

@time maxcosump()
@time limunlim()
