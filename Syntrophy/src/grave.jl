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
