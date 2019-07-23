using Syntrophy
using DifferentialEquations
# This is a script to write my testing code into
# Anything reusable should be moved into Syntrophy.jl as a seperate callable function

function npops(du,u,p::Array{Float64,1},nuts::Array{Nut,1},reacs::Array{React,1},mics::Array{Microbe,1},t::Float64)
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
    # MUST CALCULATE q

    # Do nutrients first
    for i = 1:N
        du[i] = α[i]-(δ[i]+sum(q[i,1:end]))*u[i]
    end
    # Then calculate population changes
    for i = N+1:N+M
        if E[i-N] >= 0.0 # find if growing or decaying
            # NEED TO SORT OUT DILUTION RATE HERE
            du[i] = E[i-N]*Y*u[i] - δ*u[i]
        else
            du[i] = -E[i-N]*Γ*u[i] - δ*u[i]
        end
    end
    return(du)
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
    reac = [React(1,[1,2],[1,6],[3,4],[6,6],ΔG0)]
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
    qm = 0.00 # NEED TO FIND THIS ONE
    ΔGATP = 75000 # Gibbs free energy of formation of ATP in a standard cell
    p = [Y,Γ,K,qm,ΔGATP]
    u0 = [concs;pops]
    tspan = (0.0,100.0)
    # NEED TO WORK OUT HOW TO PASS THIS PROPERLY
    prob = ODEProblem(npops,u0,p,nuts,reac,mics,tspan)
    sol = solve(prob)
    return(nothing)
end

@time main()
