# Graveyard for old functions that I should delete eventually

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

# Maybe just use a function to find threshold θ instead if we're just intrested in showing the tradeoff
function θineq(η::Float64,qm::Float64,m::Float64,kr::Float64)
    x = 0.1 # Starting with 10% as a rough guess
    θi = x*(η*qm - m)/(η*qm + kr*m)
    return(θi)
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

# Combination of the two tradeoff's to see what happens
function combtrade(set::Array{Float64,2})
    # Find length of parameter set
    L = size(set,1)
    # check that the initial parameter set is defined
    if set[1,1] <= 0.0 || set[1,2] <= 0.0 || set[1,3] <= 0.0 || set[1,4] <= 0.0
        println("Error: Initial parameter set provided incorrectly.")
    end
    # Now fill in remainder of parameter set
    for i = 2:L
        set[i,1] = set[1,1]/i
        set[i,2] = i*set[1,2]
        set[i,3] = set[1,3]/i
        set[i,4] = i*set[1,4]
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
    ps = zeros(8,20,4)
    # Asign these parameters as first parameters for each tradeoff
    for i = 1:size(ps,1)
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
    ps[8,:,:] = combtrade(ps[8,:,:])
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
    plot(r[8,:]*1.0e17,N[8,:])
    savefig("Output/CombrvsK.png")
    plot(r[8,:]*1.0e17,θt[8,:])
    savefig("Output/Combrvsθ.png")
    plot(N[8,:],θt[8,:])
    savefig("Output/CombKvsθ.png")
    return(nothing)
end

# This is a nice function that took me a bit of time to write so I want to keep
# function to create alternate parameters from an initial parameter set stored as first row
function shiftks(kset::Array{Float64,2},N::Int64,h::Float64,maxr::Float64)
    # Check if array is appropriate size
    if size(kset,1) != N
        println("Error: Parameter array is too small. Expected $(N) rows and got $(size(kset,1)).")
        error()
    end
    # Generates parameters set at various angles
    for i = 2:N
        # One parameter
        c1 = cos(2*(i-1)*pi/(N-1))
        mc1 = abs(c1)
        dc1 = sign(c1)
        if dc1 == 1.0
            # The 1+mc1*(h-1) keeps the scaling sensible
            kset[i,1] = kset[1,1]*(1+mc1*(h-1))
        else
            kset[i,1] = kset[1,1]/(1+mc1*(h-1))
        end
        # And then the other
        c2 = sin(2*(i-1)*pi/(N-1))
        mc2 = abs(c2)
        dc2 = sign(c2)
        if dc2 == 1.0
            kset[i,3] = kset[1,3]*(1+mc2*(h-1))
        else
            kset[i,3] = kset[1,3]/(1+mc2*(h-1))
        end
        # overwrite rates if greater than maxrate
        if kset[i,1] > maxr
            kset[i,1] = maxr
        end
        if kset[i,3] > maxr
            kset[i,3] = maxr
        end
    end
    return(kset)
end

# Function to find K2 from three other rates plus the equilbrium constant, plus checks if reasonable
function findK2(k1::Float64,k2::Float64,K1::Float64,KeQ::Float64,maxr::Float64)
    K2 = (k1*k2)/(K1*KeQ)
    # Test if K2 is a reasonable rate and if not return Inf
    if K2 >= maxr
        return(Inf)
    else
        return(K2)
    end
    return(K2)
end

# function to create alternate parameters from an initial parameter set stored as first row
function shiftks(kset::Array{Float64,2},h::Float64)
    # Don't change first line, just use to update the second and third lines
    kset[2,3] = h*kset[1,3]
    kset[3,3] = kset[1,3]/h
    return(kset)
end

# function to find the value of KS that maximises
function maxK(k1::Float64,k2::Float64,KeQ::Float64,maxr::Float64)
    # Starting values of rates
    K1 = 1.00e3
    # factor to increase/decrease by
    h = 1.001
    # Make array to store parameters
    N = 3 # testing 8+1 parameter sets
    kset = fill(Inf,(N,4))
    kset[:,2] .= k2
    kset[:,1] .= k1
    kset[1,3] = K1
    KS = fill(Inf,N)
    # Now start while loop
    fsh = false
    while fsh == false
        # Use first row of parameter set to generate other rows
        kset = shiftks(kset,h)
        # Find and add K2 values
        for i = 1:N
            kset[i,4] = findK2(kset[i,1],kset[i,2],kset[i,3],KeQ,maxr)
        end
        # Calculate KS values for the sets
        for i = 1:N
            if isfinite(kset[i,4])
                # Then calculate values for KS
                KS[i] = satK(kset[i,1],kset[i,2],kset[i,3])
            else
                # If K2 value is set Inf then ignore set
                KS[i] = Inf
            end
        end
        # Choose smallest value and make this the new parameter set
        I = argmin(KS)
        # Now overwrite
        kset[1,:] .= kset[I,:]
        # Stop while loop if best parameter set is self
        if I == 1
            fsh = true
        end
    end
    K1 = kset[1,3]
    K2 = kset[1,4]
    return(K1,K2)
end

# Function to do some vague investigating into the tradeoffs
function tradeinvest()
    # First need to define some basic parameters
    # Bunch of assumptions about the environment
    α = 5.55e-6 # These two rates can be changed
    δ = 2.00e-4
    θ = 0.0 # Thermodynamic inhibition not reached
    CO = 0.21 # High value so oxegen isn't limiting
    # Will consider glucose to begin with and then move onto lower free energy changes
    ΔG0 = -2843800.0
    # Assume 38 ATP generated per glucose metabolised
    η = 38.0
    m = 2.16e-19 # maintainance, seems fine shouldn't alter results too much
    Y = 2.36e13 # yield in cells per mole of ATP, scaling constant might need to be redone, but shouldn't alter results
    E0 = 2.5e-20 # Moles of enzyme per cell, not 100% sure of this number, should be improved later
    # Define other thermodynamically relevant constants
    ΔGATP = 75000.0 # Gibbs free energy of formation of ATP in a standard cell
    Temp = 312.0 # Temperature that growth is occuring at in Kelvin
    # Use to find equilbrium constant
    KeQ = Keq(ΔG0,η,ΔGATP,Temp)
    # Now want to choose maximum populations for a range of k values
    k2s = collect(10.0:10.0:1000.0)
    qm = zeros(length(k2s))
    KS = zeros(length(k2s))
    kr = zeros(length(k2s))
    QT = zeros(length(k2s))
    Ns = zeros(length(k2s))
    maxr = 1.00e6
    k1 = 1.00e4 # Now a fixed quantity, should be carefully chosen
    for i = 1:length(k2s)
        k2 = k2s[i]
        K1, K2 = maxK(k1,k2,KeQ,maxr)
        qm[i], KS[i], _, kr[i] = qKmY(k1,K1,k2,K2,E0)
        QT[i] = Qineq(η,qm[i],m,kr[i],KeQ)
        _, _, Ns[i] = stead(KS[i],kr[i],η,qm[i],m,CO,α,δ,θ)
    end
    pyplot(dpi=200)
    Lqm = L"q_m"
    plot(qm*1e18,KS,xlabel="$(Lqm) (10^-18)")
    savefig("Output/qmvsKS.png")
    plot(qm*1e18,kr,xlabel="$(Lqm) (10^-18)")
    savefig("Output/qmvskr.png")
    plot(qm*1e18,Ns,xlabel="$(Lqm) (10^-18)")
    savefig("Output/qmvsNs.png")
    plot(qm*1e18,QT,xlabel="$(Lqm) (10^-18)")
    savefig("Output/qmvsQT.png")
    return(nothing)
end
