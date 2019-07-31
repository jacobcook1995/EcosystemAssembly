using Syntrophy
using DifferentialEquations
using Plots

# This is a script to write my testing code into
# Anything reusable should be moved into Syntrophy.jl as a seperate callable function

# Function to update population and nutrient concentrations
# Throughout this function there are sections that are very specific needs to be further generalised!
function npops(du::Array{Float64,1},u::Array{Float64,1},p::Array{Float64,1},nuts::Array{Nut,1},reacs::Array{React,1},
                mics::Array{Microbe,1},ex::Array{Int64,1},t::Float64)
    # Extract required parameters
    Y = p[1]
    Γ = p[2]
    # Extract relevant data from vector of nutrients
    α = nuts.↦:α
    δ = nuts.↦:δ
    N = length(nuts) # Number of nutrients
    # And relevant data from vector of microbes
    η = mics.↦:η
    m = mics.↦:m
    δ2 = mics.↦:δ
    M = length(mics) # Number of microbes
    # Extract reaction stochiometry
    stc = (reacs.↦:stc)[1]
    ΔG0 = (reacs.↦:ΔG0)[1]
    # Now calculate q
    q = thermrate(u[1:N],u[N+1:end],p[3],p[4],p[5],ΔG0,p[6],stc,η)
    # Make vector to store nutrient changes due to consumption
    δX = zeros(N)
    for i = 1:length(stc)
        for j = N+1:N+M
            δX[i] += stc[i]*q[j-N]*u[j] # Assumes single reaction
        end
    end
    # q has no dependance on population
    # Now update nutrients
    for i = 1:N
        if (i in ex) == false
            du[i] = α[i]-δ[i]*u[i]+δX[i]
        else
            du[i] = 0
        end
    end
    # Then calculate population changes
    for i = N+1:N+M
        j = i-N
        E = netE(η[j],q[j],m[j])
        if E >= 0.0 # find if growing or decaying
            du[i] = E*Y*u[i] - δ2[j]*u[i]
        else
            du[i] = E*Γ*u[i] - δ2[j]*u[i]
        end
    end
    return(du)
end

# function to calculate the thermodynamic consumption rate
# This function also needs genralising
function thermrate(concs::Array{Float64,1},pops::Array{Float64,1},K::Float64,qm::Float64,ΔGATP::Float64,
                    ΔG0::Float64,Temp::Float64,stoc::Array{Int64,1},η::Array{Float64,1})
    # concs => Vector of nutrient concentrations
    # pops => Vector of population densities
    # K => Saturation constant for the substrate
    # qm => Maximal reaction rate for substrate
    # stoc => stochiometry vector
    # ΔGATP => Gibbs free energy to form ATP in standard cell
    # ΔG0 => Standard gibbs free energy of reaction
    # Temp => Temperature in Kelvin
    # η => free energy use strategy, mol of ATP per mol of substrate
    ############ START OF FUNCTION ###################

    # Initialise consumption matrix
    q = zeros(length(pops))

    # Calculate substrate coefficent
    S = 1
    for i = 1:length(stoc)
        if stoc[i] < 0
            S *= concs[i]^(-stoc[i])
        end
    end
    # Loop over population to find rates q
    for i = 1:length(pops)
        # Call function to find thermodynamic factor θ
        θ = θT(concs,stoc,ΔGATP,ΔG0,η[i],Temp)
        # Only η changes between species
        q[i] = qm*S*(1-θ)/(K+S*(1+θ))
    end
    return(q)
end

function gluc()
    # Nutrient variables
    α = 3.00*10^(-8)
    δ = 1.00*10^(-6) # Death rate that doesn't wash out cells
    # make nutrients
    # 1 = glucose, 2 = oxegen, 3 = bicarbonate, 4 = hydrogen ion
    nuts = [Nut(1,α,δ),Nut(2,6*α,δ),Nut(3,0,δ),Nut(4,0,δ)]
    # Now make reactions
    ΔG0 = -2843800.0
    reac = [React(1,[1,2,3,4],[-1,-6,6,6],ΔG0)]
    # microbe variables
    η1 = 41.2
    η2 = 41.0
    # maintainance equal
    m = 2.16*10^(-19)
    # both have same substrate and end product
    r = 1
    # make microbes
    mics = [Microbe(η1,m,r,δ),Microbe(η2,m,r,δ)]
    # Set intial populations and nutrient concentrations
    pops = 10.0*ones(length(mics))
    concs = zeros(length(nuts))
    concs[1] = 0.3# start with high amount of glucose
    concs[2] = 0.0018 # WHY NOT JUST FIX O2 and pH?
    concs[3] = 1.00*10^(-9)
    concs[4] = 1.00*10^(-7) # pH 7
    # Define some constants
    Y = 2.36*10^(13) # biomass yield cells per mol ATP
    Γ = 1.16*10^(12) # starvation rate cells per mol ATP (deficit)
    K = 2.00*10^(-10) # Saturation constant mol ml^(−1)
    qm = 4.44*10^(-13) # Maximal possible growth rate mol cell^(-1) s^(-1)
    ΔGATP = 75000.0 # Gibbs free energy of formation of ATP in a standard cell
    Temp = 312.0 # Temperature that growth is occuring at in Kelvin
    p = [Y,Γ,K,qm,ΔGATP,Temp]
    u0 = [concs;pops]
    tspan = (0.0,2500000.0)
    ex = [2,4]
    # Make reduced version of function inputting unchanging microbes
    f(du,u,p,t) = npops(du,u,p,nuts,reac,mics,ex,t)
    prob = ODEProblem(f,u0,tspan,p)
    sol = solve(prob,adaptive=false,dt=2000) # turned dt down to make plots look nicer
    # Now do plotting
    plot(sol.t,sol'[:,5:6],label=["\\eta = $((mics.↦:η)[1])","\\eta = $((mics.↦:η)[2])"])
    savefig("Output/Populations$(η1)vs$(η2).png")
    return(nothing)
end

@time gluc()
