using Syntrophy
using Plots
using DifferentialEquations

# This is a script to write my testing code into
# Anything reusable should be moved into Syntrophy.jl as a seperate callable function

# Function to update population and nutrient concentrations
# This is run for a single population utilising a single reaction
function singlepop(du::Array{Float64,1},u::Array{Float64,1},p::Array{Float64,1},nuts::Array{Nut,1},reacs::Array{React,1},
                mics::Array{Microbe,1},t::Float64)
    # Extract required parameters
    Y = p[1]
    # Extract relevant data from vector of nutrients
    α = nuts.↦:α
    δ = nuts.↦:δ
    con = nuts.↦:cst
    N = length(nuts) # Number of nutrients
    # And relevant data from vector of microbes
    η = (mics.↦:η)[1]
    m = (mics.↦:m)[1] # running for single microbe
    M = length(mics) # Number of microbes
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

# function to calculate the rate of substrate consumption q
function qrate(concs::Array{Float64,1},KS::Float64,qm::Float64,ΔGATP::Float64,
                    ΔG0::Float64,Temp::Float64,stoc::Array{Int64,1},η::Float64)
    # concs => Vector of nutrient concentrations
    # KS => Saturation constant for the substrate
    # qm => Maximal reaction rate for substrate
    # stoc => stochiometry vector
    # ΔGATP => Gibbs free energy to form ATP in standard cell
    # ΔG0 => Standard gibbs free energy of reaction
    # Temp => Temperature in Kelvin
    # η => free energy use strategy, mol of ATP per mol of substrate
    ############ START OF FUNCTION ###################

    # calulate substrate coefficent
    S = SCoef(concs,stoc)
    # Call function to find thermodynamic factor θ
    θ = θT(concs,stoc,ΔGATP,ΔG0,η,Temp)
    # Only η changes between species
    q = qm*S*(1-θ)/(KS+S*(1+θ))
    # Catch unbiological negative rate case
    if q < 0.0
        q = 0.0
    end
    return(q)
end



function gluc()
    # Nutrient variables
    α = 5.55*10^(-6)
    δ = 1.00*10^(-4)
    # make nutrients
    # 1 = glucose, 2 = oxegen, 3 = bicarbonate, 4 = hydrogen ion
    nuts = [Nut(1,false,α,δ),Nut(2,true,0,0),Nut(3,false,0,δ),Nut(4,true,0,0)]
    # Now make reactions
    ΔG0 = -2843800.0
    reac = [React(1,[1,2,3,4],[-1,-6,6,6],ΔG0)]
    # microbe variables
    η = 42.0
    r = 1 # Only reaction
    m = 2.16*10^(-19) # maintainance
    # Considering 1 microbe with no maintaince and no dilution
    mics = [Microbe(η,m,1,0.0)]
    # Set intial populations and nutrient concentrations
    pops = 100.0*ones(length(mics))
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
    θf = θT(sol'[end,1:4],stoc,ΔGATP,ΔG0,η,Temp)
    println(θf)
    return(nothing)
end

@time gluc()
