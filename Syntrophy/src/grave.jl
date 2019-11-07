# Graveyard for old functions that I should delete eventually

# Function to predict steady state populations in non-thermodynamically limited case
function predict()
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
    η = 41.0 # this is the actual physiological value
    ΔGATP = 75000.0 # Gibbs free energy of formation of ATP in a standard cell
    r = 1 # Only reaction
    # Define kinetic parameters explicitly
    E0 = 2.5*10.0^(-20) # Somewhat fudged should be right order of magnitude
    # These can be used as reference values corresponding to the case of E0ref
    m = 2.16*10.0^(-19) # maintainance
    Y = 2.36*10.0^(13) # yield in cells per mole of ATP
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
    k2 = 140.0
    k1, k2, K1, K2 = parak(k2,ΔG0,η,ΔGATP,Temp)
    println("k1 = $(k1)")
    println("k2 = $(k2)")
    println("K1 = $(K1)")
    println("K2 = $(K2)")
    # Find KS, qm, maintainance and yield using function
    qm, KS, KP, kr = qKmY(k1,K1,k2,K2,E0)
    println("qm = $(qm)")
    println("KS = $(KS)")
    println("KP = $(KP)")
    println("kr = $(kr)")
    # Considering 1 microbe with maintaince but no dilution
    mics = Microbe(η,m,1,0.0)
    p = [Y,KS,qm,ΔGATP,Temp,kr]
    # Make reduced version of function inputting unchanging microbes
    f(du,u,p,t) = singlepop(du,u,p,nuts,reac,mics,t)
    prob = ODEProblem(f,u0,tspan,p)
    sol = DifferentialEquations.solve(prob,adaptive=false,dt=500) # turned dt down to make plots look nicer
    println("K = $(maximum(sol'[:,5]))")
    println("min S = $(minimum(sol'[:,1]))") # Declines from initially high value
    println("max P = $(maximum(sol'[:,3]))")
    # Find steady state value of θ
    stoc = (reac.↦:stc)[1]
    θ = θT(sol'[end,1:4],stoc,ΔGATP,ΔG0,η,Temp)
    println("Steady state θ = $(θ)")
    # Now find and print predicted steady state values
    S, P, X = stead(KS,kr,η,qm,m,concs[2],α,δ,θ)
    println("predicted K = $(X)")
    println("predicted S = $(S)")
    println("predicted P = $(P)")
    # Same but using my more complex function
    KeQ = Keq(ΔG0,η,ΔGATP,Temp)
    # pyplot()
    # plot(sol'[:,5])
    # hline!([X])
    # savefig("Output/test1.png")
    # plot(sol'[:,3])
    # hline!([P])
    # savefig("Output/test2.png")
    # plot(sol'[:,1])
    # hline!([S])
    # savefig("Output/test3.png")
    return(nothing)
end

# function to return k parameters based on a single k value
function parak(k2::Float64,K2::Float64,ΔG0::Float64,η::Float64,ΔGATP::Float64,Temp::Float64)
    K1 = 140.0 # K1 must match value of k2 used for k_R = 1
    # Now work out equlibrium constant K in order to find final rate
    K = Keq(ΔG0,η,ΔGATP,Temp)
    k1 = K*K1*K2/(k2)
    return(k1,k2,K1,K2)
end

# function to return k parameters based on a single k value, now for thermodynamic tradeoff
function parakT(k2::Float64,k1::Float64,ΔG0::Float64,η::Float64,ΔGATP::Float64,Temp::Float64)
    K1 = 140.0 # K1 must match value of k2 used for k_R = 1
    # Now work out equlibrium constant K in order to find final rate
    K = Keq(ΔG0,η,ΔGATP,Temp)
    K2 = k1*k2/(K*K1)
    return(k1,k2,K1,K2)
end

# function to find K2 based on changed value of η
function newK2(ΔG0::Float64,η::Float64,ΔGATP::Float64,Temp::Float64)
    K2 = 1.28*10.0^(8) # Can ignore other rates as base case has k_R = 1
    # Now work out equlibrium constant K in for η = 38
    Kb = Keq(ΔG0,38.0,ΔGATP,Temp)
    # Now find equlbrium constant for η we are considering
    Ka = Keq(ΔG0,η,ΔGATP,Temp)
    K2 *= (Kb/Ka)
    return(K2)
end

# Calculate maxrate for direct tradeoff
function maxrate(k2::Float64,ΔG0::Float64,η::Float64,ΔGATP::Float64,Temp::Float64,E0::Float64,Y::Float64,f::Function,u0::Array{Float64,1})
    if η != 38.0
        # First need find new value of K2
        K2 = newK2(ΔG0,η,ΔGATP,Temp)
    else
        K2 = 1.28*10.0^(8) # Fix value of K_2
    end
    # Calculate other rates to match k
    k1, k2, K1, K2 = parak(k2,K2,ΔG0,η,ΔGATP,Temp)
    # Use rates to obtain required parameters
    qm, KS, KP, kr = qKmY(k1,K1,k2,K2,E0)
    # Now need to calculate maximal rate mr
    p = [Y,KS,qm,ΔGATP,Temp,kr] # collect parameters
    # put parameters into function
    tspan = (0.0,5000000.0)
    prob = ODEProblem(f,u0,tspan,p)
    # then solve
    sol = solve(prob,adaptive=false,dt=100) # Very detailed
    # Need a check that population has ceased growing
    # Check if pop has changed more than 0.1% in last 10 time steps
    diff = (sol'[end,5]-sol'[end-10,5])/sol'[end,5]
    if diff > 0.001
        println("Not a stable population!")
    end
    # Find final (maximum) population
    mp = sol'[end,5]
    # Then need to find actual q rate along trajectory
    qs = zeros(length(sol.t))
    for i = 1:length(qs)
        qs[i] = qrate(sol'[i,1:4],KS,qm,ΔGATP,ΔG0,Temp,[-1,-6,6,6],η,kr)
    end
    # Take max of this new vector
    mr = maximum(qs)
    # return both qm and actual maximal observed rate
    return(qm,mr,mp,KS)
end

# Calculate maxrate for indirect thermodynamic (T) tradeoff
function maxrateT(k2::Float64,ΔG0::Float64,η::Float64,ΔGATP::Float64,Temp::Float64,E0::Float64,Y::Float64,f::Function,u0::Array{Float64,1})
    k1 = 1.17*10.0^(7) # Fix value of k_1
    # Calculate other rates to match k
    k1, k2, K1, K2 = parakT(k2,k1,ΔG0,η,ΔGATP,Temp)
    # Use rates to obtain required parameters
    qm, KS, KP, kr = qKmY(k1,K1,k2,K2,E0)
    # Now need to calculate maximal rate mr
    p = [Y,KS,qm,ΔGATP,Temp,kr] # collect parameters
    # put parameters into function
    tspan = (0.0,5000000.0)
    prob = ODEProblem(f,u0,tspan,p)
    # then solve
    sol = solve(prob,adaptive=false,dt=100) # Very detailed
    # Need a check that population has ceased growing
    # Check if pop has changed more than 0.1% in last 10 time steps
    diff = (sol'[end,5]-sol'[end-10,5])/sol'[end,5]
    if diff > 0.001
        println("Not a stable population!")
    end
    # Find final (maximum) population
    mp = sol'[end,5]
    # Then need to find actual q rate along trajectory
    qs = zeros(length(sol.t))
    for i = 1:length(qs)
        qs[i] = qrate(sol'[i,1:4],KS,qm,ΔGATP,ΔG0,Temp,[-1,-6,6,6],η,kr)
    end
    # Take max of this new vector
    mr = maximum(qs)
    # Find reaction quotient
    Q = QCoef(sol'[end,1:4],[-1,-6,6,6])
    # return both qm and actual maximal observed rate
    return(qm,mr,mp,KS,kr,Q)
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
    u0[1] = 5.55 # high initial concentration to ensure growth near maximum
    u0[2] = 0.21 # High value so oxegen isn't limiting
    u0[3] = 0.0 # No initial concentration
    u0[4] = 1.00*10.0^(-7) # pH 7
    u0[5] = 100.0
    # Preallocate vectors
    qm = zeros(20)
    mr = zeros(length(qm))
    mp = zeros(length(qm))
    KS = zeros(length(qm))
    # In this case where η = 38 use a preset value of K2
    # Loop over vectors
    for i = 1:length(qm)
        # We now want to increase k_{+2} to increase q_m
        k2 = 140.0*i
        qm[i], mr[i], mp[i], KS[i] = maxrate(k2,ΔG0,η,ΔGATP,Temp,E0,Y,f,u0)
    end
    # Switch backends
    pyplot(dpi=200)
    Ns = L"N^{\ast}"
    p14 = L"10^{14}"
    plot(qm*10.0^17,mp*10.0^(-14),xlabel=L"q_m\;(s^{-1}\,10^{-17})",ylabel="$(Ns) (cells $(p14))")
    savefig("Output/qmvsmp.png")
    plot(qm*10.0^17,mr*10.0^19,xlabel=L"q_m\;(s^{-1}\,10^{-17})",ylabel=L"q\;(s^{-1}\,10^{-19})")
    savefig("Output/qmvsmr.png")
    return(nothing)
end

# function to test how the previous tradeoff changes in thermodynamic limit
function testq2()
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
    # thermodynamically limiting η
    η = 41.25
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
    u0[1] = 0.0555 # high initial concentration to ensure growth near maximum
    u0[2] = 0.21 # High value so oxegen isn't limiting
    u0[3] = 0.0 # No initial concentration
    u0[4] = 1.00*10.0^(-7) # pH 7
    u0[5] = 100.0
    # Preallocate vectors
    qm = zeros(20)
    mr = zeros(length(qm))
    mp = zeros(length(qm))
    KS = zeros(length(qm))
    Nst = zeros(length(qm))
    # Loop over vectors
    for i = 1:length(qm)
        # We now want to increase k_{+2} to increase q_m
        k2 = 140.0*i
        qm[i], mr[i], mp[i], KS[i] = maxrate(k2,ΔG0,η,ΔGATP,Temp,E0,Y,f,u0)
        Nst[i] = (η/m)*(α - δ*m*KS[i]/((η*qm[i] - m)*(u0[2]^6)))
    end
    # Switch backends
    pyplot(dpi=200)
    Ns = L"N^{\ast}"
    p14 = L"10^{14}"
    plot(qm*10.0^(17),mp*10.0^(-14),xlabel=L"q_m\;(s^{-1}\,10^{-17})",ylabel="$(Ns) (cells $(p14))")
    plot!(qm*10.0^(17),Nst*10.0^(-14))
    savefig("Output/Tqmvsmp$(η).png")
    plot(qm*10.0^(17),mr*10.0^19,xlabel=L"q_m\;(s^{-1}\,10^{-17})",ylabel=L"q\;(s^{-1}\,10^{-19})")
    savefig("Output/Tqmvsmr.png")
    return(nothing)
end

# function to test other tradeoff away from thermodynamic limitation
function testq3()
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
    # physiological value of η = 38.0, change to investigate thermodynamic inhibition
    η = 41.1 # 38.0
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
    u0[1] = 0.0555 # high initial concentration to ensure growth near maximum
    u0[2] = 0.21 # High value so oxegen isn't limiting
    u0[3] = 0.0 # No initial concentration
    u0[4] = 1.00*10.0^(-7) # pH 7
    u0[5] = 100.0
    # Preallocate vectors
    qm = zeros(20)
    mr = zeros(length(qm))
    mp = zeros(length(qm))
    KS = zeros(length(qm))
    kr = zeros(length(qm))
    Nst = zeros(length(qm))
    Qp = zeros(length(qm))
    Qa = zeros(length(qm))
    # Loop over vectors
    for i = 1:length(qm)
        println("Step $(i):")
        # We now want to increase k_{+2} to increase q_m
        k2 = 140.0*i
        qm[i], mr[i], mp[i], KS[i], kr[i], Qa[i] = maxrateT(k2,ΔG0,η,ΔGATP,Temp,E0,Y,f,u0)
        Nst[i] = (η/m)*(α - δ*m*KS[i]/((η*qm[i] - m)*(u0[2]^6)))
        KeQ = Keq(ΔG0,η,ΔGATP,Temp)
        Qp[i] = Qineq(η,qm[i],m,kr[i],KeQ)
    end
    # Switch backends
    pyplot(dpi=200)
    Ns = L"N^{\ast}"
    p14 = L"10^{14}"
    plot(qm*10.0^(17),mp*10.0^(-14),xlabel=L"q_m\;(s^{-1}\,10^{-17})",ylabel="$(Ns) (cells $(p14))")
    plot!(qm*10.0^(17),Nst*10.0^(-14))
    savefig("Output/RatevsPopT$(η).png")
    println("Difference vs predicted:")
    # println((Nst.-mp)*10.0^(-14))
    # println((Nst.-mp)./Nst)
    plot(qm*10.0^(17),mr*10.0^19,xlabel=L"q_m\;(s^{-1}\,10^{-17})",ylabel=L"q\;(s^{-1}\,10^{-19})")
    savefig("Output/ActualvsPredictedRateT$(η).png")
    plot(qm*10.0^(17),Qp*10.0^41,label="predict")
    # plot!(qm*10.0^(17),Qa*10.0^41,label="actual")
    savefig("Output/ActualvsPredictedQ$(η).png")
    return(nothing)
end


# Quick test function to test my other tradeoffs
function test()
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
    # physiological value of η = 38.0, change to investigate thermodynamic inhibition
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
    # Set constant chemical species
    u0 = zeros(2)
    u0[1] = 0.21 # High value so oxegen isn't limiting
    u0[2] = 1.00*10.0^(-7) # pH 7
    # Define parameters by hand
    fac = 1#1000000
    K1 = 140.0*fac
    k2 = 140.0*fac
    k1 = 1.17e7
    K2 = 1.28e8
    # Use to calculate parameters
    qm, KS, KP, kr = qKmY(k1,K1,k2,K2,E0)
    println("qm = $(qm)")
    println("KS = $(KS)")
    println("KP = $(KP)")
    println("kr = $(kr)")
    # Use to find steady state assuming θ≈0
    S, P, X = stead(KS,kr,η,qm,m,u0[1],α,δ,0.0)
    # Then calculate actula θ of this steady state
    stoc = (reac.↦:stc)[1]
    θt = θT([S,u0[1],P,u0[2]],stoc,ΔGATP,ΔG0,η,Temp)
end
