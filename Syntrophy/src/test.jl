using Syntrophy
using DifferentialEquations
# This is a script to write my testing code into
# Anything reusable should be moved into Syntrophy.jl as a seperate callable function

# Function to update population and nutrient concentrations
# Throughout this function there are sections that are very specific needs to be further generalised!
function npops(du::Array{Float64,1},u::Array{Float64,1},p::Array{Float64,1},nuts::Array{Nut,1},reacs::Array{React,1},mics::Array{Microbe,1},t::Float64)
    # Extract required parameters
    Y = p[1]
    Γ = p[2]
    # Extract relevant data from vector of nutrients
    α = nuts.↦:α
    δ = nuts.↦:δ
    N = length(nuts) # Number of nutrients
    # And relevant data from vector of microbes
    η = mics.↦:η
    M = length(mics) # Number of microbes
    # Now calculate q
    q = thermrate(u[1:N],u[N+1:end],p[3],p[4])
    # Make vector to store nutrient changes due to consumption
    δX = zeros(N)
    # Extract reaction stochiometry
    stc = (reacs.↦:stc)[1]
    for i = 1:length(stc)
        δX[i] = stc[i]*sum(q) # Assumes single reaction
    end
    println(q)
    println(δX)
    error()
    # Do nutrients first
    for i = 1:N
        du[i] = α[i]-(δ[i]+δX[i])*u[i]
    end
    # Then calculate population changes
    for i = N+1:N+M
        if E[i-N] >= 0.0 # find if growing or decaying
            # NEED TO SORT OUT DILUTION RATE HERE
            du[i] = E[i-N]*Y*u[i] - δ[i]*u[i]
        else
            du[i] = -E[i-N]*Γ*u[i] - δ[i]*u[i]
        end
    end
    return(du)
end

# function to calculate the thermodynamic consumption rate
# This function also needs genralising
function thermrate(concs::Array{Float64,1},pops::Array{Float64,1},K::Float64,qm::Float64,stoc::Array{Int64,1})
    # concs => Vector of nutrient concentrations
    # pops => Vector of population densities
    # K => Saturation constant for the substrate
    # qm => Maximal reaction rate for substrate
    # stoc => stochiometry vector
    ############ START OF FUNCTION ###################

    # Initialise consumption matrix
    q = zeros(length(pops))

    # Calculate substrate coefficent
    S = 0
    for i = 1:length(stoc)
        if stoc[i] < 0
            S *= concs[i]^(-stoc[i])
        end
    end

    # Loop over population to find rates q
    for i = 1:length(pops)
        # Call function to find thermodynamic factor θ
        θ = θT(concs,stoc,ΔGATP,ΔG0,η,Temp)
        # NEED TO ADD IN ΔGATP,ΔG0,η,Temp
        # ONLY η changes between species
        q[i] = qm*S*(1-θ[i])/(K+S*(1+θ[i]))
    end
    return(q)
end

function main()
    # Nutrient variables
    α = 4.70*10^(-9)
    δ = 3.50*10^(-5)
    # make nutrients
    # 1 = glucose, 2 = oxegen, 3 = bicarbonate, 4 = hydrogen ion
    nuts = [Nut(1,α,δ),Nut(2,6*α,δ),Nut(3,0,δ),Nut(4,0,0)]
    # Now make reactions
    ΔG0 = -2843800.0
    reac = [React(1,[1,2,3,4],[-1,-6,6,6],ΔG0)]
    # microbe variables
    η1 = 1.0
    η2 = 0.5
    # maintainance equal
    m = 2.16*10^(-19)
    # both have same substrate and end product
    r = 1
    # make microbes
    mics = [Microbe(η1,m,r),Microbe(η2,m,r)]
    # Set intial populations and nutrient concentrations
    pops = 0.001*ones(length(mics))
    concs = zeros(length(nuts))
    concs[1] = 1.00*10^(-9)
    concs[2] = 2.70*10^(-7)
    concs[3] = 1.00*10^(-9)
    concs[4] = 1.00*10^(-7) # pH 7
    # Define some constants
    Y = 2.36*10^(13) # biomass yield cells per mol ATP
    Γ = 1.16*10^(12) # starvation rate cells per mol ATP (deficit)
    K = 2.00*10^(-10) # Saturation constant mol ml^(−1)
    qm = 4.44*10^(-13) # Maximal possible growth rate mol cell^(-1) s^(-1)
    ΔGATP = 75000 # Gibbs free energy of formation of ATP in a standard cell
    p = [Y,Γ,K,qm,ΔGATP]
    u0 = [concs;pops]
    tspan = (0.0,100.0)
    # Make reduced version of function inputting unchanging microbes
    f(du,u,p,t) = npops(du,u,p,nuts,reac,mics,t)
    prob = ODEProblem(f,u0,tspan,p)
    sol = solve(prob)
    return(nothing)
end

@time main()
