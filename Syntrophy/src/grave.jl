# Graveyard for old functions that I should delete eventually

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
