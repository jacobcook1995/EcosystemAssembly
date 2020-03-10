# Script that makes the parameters needed for different types of simulation
# of the model with inhibition
using Assembly

export initialise, initialise_η, initialise_η2

# function to generate a set of random reactions for the model
function rand_reactions(O::Int64,M::Int64,μrange::Float64,T::Float64)
    # make uniform distribution d
    h = 6.0
    l = 1.0
    d = Uniform(l,h)
    # preallocate output
    RP = zeros(Int64,O,2)
    ΔG = zeros(O)
    # Preallocate vector to store chemical potential at standard conditions
    μ0 = zeros(M)
    # Find random intervals
    rs = rand(M-1)
    rs = rs/sum(rs)
    # Take maximum value to be μrange, this is unphysical but works as we are intrested in changes
    μ0[1] = μrange
    for i = 2:M
        μ0[i] = μ0[i-1] - rs[i-1]*μrange
    end
    # Now randomly assign reactions
    i = 0
    while i < O
        # Randomly choose reactant
        r = rand(1:M)
        # Randomly choose product
        p = rand(1:M)
        # Check that r and p are not already stored
        if all((r .∉ RP[:,1]) .| (p .∉ RP[:,2]))
            # find Gibbs free energy change
            dG = μ0[p] - μ0[r]
            # Find equlibrium constant when 1 ATP is generated
            K = Keq(T,1.0,dG)
            # Find log base 10 of this equilbrium constant
            logK = log(10,K)
            # Accept if equilbrium ratio is less than 10:1 in favour of reactant
            if logK > -l
                i += 1
                RP[i,1] = r
                RP[i,2] = p
                ΔG[i] = dG
            # Disregard if too heavily in favour of reactant
            elseif logK < -h
            # Otherwise randomly decide whether to accept
            else
                # draw random number between 1.0 and 6.0
                rn = rand(d)
                # Accept if this is
                if rn > -logK
                    i += 1
                    RP[i,1] = r
                    RP[i,2] = p
                    ΔG[i] = dG
                end
            end
        end
    end
    return(RP,ΔG)
end

# function to randomly choose the reactions that a specific microbe can make use of
function choose_reactions(O::Int64,mR::Float64,sdR::Float64)
    @assert mR - 5*sdR >= 0.0 "This choice could result in a negative number of reactions per microbe"
    # Make required Gaussian distribution using the provided mean (mm) and SD (sdm)
    d = Normal(mR,sdR)
    # randomly choose the number
    R = round(Int64,rand(d))
    # Preallocate vector of reaction identities
    Reacs = zeros(Int64,R)
    # Choose random first value
    Reacs[1] = rand(1:O)
    # Then fill out later values
    for i = 2:R
        good = false
        while good == false
            r = rand(1:O)
            # Check to avoid repeated values
            if r ∉ Reacs[1:i-1]
                Reacs[i] = r
                good = true
            end
        end
    end
    return(R,Reacs)
end

# Function to choose kinetic parameters, do by guassian perturbation to mean at the moment
function choose_kinetic(R::Int64,mq::Float64,sdq::Float64,mK::Float64,sdK::Float64,mk::Float64,sdk::Float64)
    @assert mq - 5*sdq >= 0.0 "This choice could result in a negative maximal rate"
    @assert mK - 5*sdK >= 0.0 "This choice could result in a negative saturation constant"
    @assert mk - 5*sdk >= 0.0 "This choice could result in a negative reversibility"
    # preallocate vectors to store kinetic paarmeters
    qm = zeros(R)
    KS = zeros(R)
    kr = zeros(R)
    # Set up distributions
    d1 = Normal(mq,sdq)
    d2 = Normal(mK,sdK)
    d3 = Normal(mk,sdk)
    # Now loop over all the reactions
    for i = 1:R
        qm[i] = rand(d1)
        KS[i] = rand(d2)
        kr[i] = rand(d3)
    end
    return(qm,KS,kr)
end

# function to choose η values for each reaction based on Gibbs free energy changes
function choose_ηs(reacs::Array{Reaction,1},Reacs::Array{Int64,1},T::Float64)
    # Preallocate memory to store η's
    η = zeros(length(Reacs))
    # Set a constant lower bound
    ηl = 1/3
    # Set minimum equilibirum product to substrate ratio
    mratio = 1e-5
    # Make beta distribution for later, parameters chosen so that distribution skews right
    d = Beta(5,1)
    for i = 1:length(η)
        # Indentify which reaction we are considering
        I = Reacs[i]
        # Find corresponding Gibbs free energy change
        dG = reacs[I].ΔG0
        # And use to determine an upper bound on η
        ηh = -(dG + Rgas*T*log(mratio))/(ΔGATP)
        η[i] = (ηh-ηl)*rand(d) + ηl
    end
    return(η)
end

# function to generate parameter set for the model with inhibition
function initialise(N::Int64,M::Int64,O::Int64,mR::Float64,sdR::Float64,mq::Float64,sdq::Float64,mK::Float64,sdK::Float64,mk::Float64,sdk::Float64)
    @assert O > mR + 5*sdR "Not enough reactions to ensure that microbes have on average mR reactions"
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # And all proportionality constants are the same for simplicity
    g = ones(N)
    # All but resource 1 is not supplied
    κ = zeros(M)
    κ[1] = 100.0
    # Assume that all δ's are equal
    δi = 1.0
    δ = δi*ones(M)
    # Find m using a function that gives a Guassian offset
    mm = 1.0
    sdm = 0.1
    m = mvector(N,mm,sdm)
    # Generate random set of reactions
    μrange = 3e6 # Chosen to be slightly larger than glucose respiration
    RP, ΔG = rand_reactions(O,M,μrange,T)
    # Preallocate vector of reactions
    reacs = Array{Reaction,1}(undef,O)
    for i = 1:O
        reacs[i] = make_Reaction(i,RP[i,1],RP[i,2],ΔG[i])
    end
    # Preallocate vector of microbes
    mics = Array{Microbe,1}(undef,N)
    # Then construct microbes
    for i = 1:N
        # For each microbe generate random set of reactions
        R, Reacs = choose_reactions(O,mR,sdR)
        # Find corresponding kinetic parameters for these reactions
        qm, KS, kr = choose_kinetic(R,mq,sdq,mK,sdK,mk,sdk)
        # Find corresponding η's for these reactions
        η = choose_ηs(reacs,Reacs,T)
        # Can finally generate microbe
        mics[i] = make_Microbe(m[i],g[i],R,Reacs,η,qm,KS,kr)
    end
    # Now make the parameter set
    ps = make_InhibParameters(N,M,O,T,κ,δ,reacs,mics)
    return(ps)
end

# function to choose single η value deterministically based on microbe number i
function single_ηs(reacs::Array{Reaction,1},T::Float64,N::Int64,i::Int64)
    # Set a constant lower bound, higher in this case as non longer sampling from right skewed distribution
    ηl = 1
    # Set minimum equilibrium product to substrate ratio
    mratio = 1e-5
    # Find corresponding Gibbs free energy change
    dG = reacs[1].ΔG0
    # And use to determine an upper bound on η
    ηh = -(dG + Rgas*T*log(mratio))/(ΔGATP)
    η = (ηh-ηl)*(i/N) + ηl
    return(η)
end

# function to generate parameter set for the case of η competition on one reaction
function initialise_η(N::Int64,mq::Float64,mK::Float64,mk::Float64)
    # Only two metabolities alterable by one reaction in this case
    M = 2
    O = 1
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # And all proportionality constants are the same for simplicity
    g = ones(N)
    # All but resource 1 is not supplied
    κ = zeros(M)
    κ[1] = 100.0
    # Assume that all δ's are equal
    δi = 1.0
    δ = δi*ones(M)
    # Use deterministic maintenance rates
    mm = 1.0
    m = mm*ones(N)
    # Generate random set of reactions
    μrange = 3e5 # Smaller μrange in this case to show η competition
    RP, ΔG = rand_reactions(O,M,μrange,T)
    # Preallocate vector of reactions
    reacs = Array{Reaction,1}(undef,O)
    for i = 1:O
        reacs[i] = make_Reaction(i,RP[i,1],RP[i,2],ΔG[i])
    end
    # Preallocate vector of microbes
    mics = Array{Microbe,1}(undef,N)
    # Then construct microbes
    for i = 1:N
        # Only one reaction which each microbe uses
        R = 1
        Reacs = [1]
        # Kinetic parameters are fixed
        qm = [mq]
        KS = [mK]
        kr = [mk]
        # Find corresponding η's for these reactions
        η = [single_ηs(reacs,T,N,i)]
        # Can finally generate microbe
        mics[i] = make_Microbe(m[i],g[i],R,Reacs,η,qm,KS,kr)
    end
    # Now make the parameter set
    ps = make_InhibParameters(N,M,O,T,κ,δ,reacs,mics)
    return(ps)
end

# function to deterministically generate 2 reactions in 3 metabolite case
function two_reactions(μrange::Float64)
    # Two reactions from 3 metabolites
    O = 2
    M = 3
    # preallocate output
    RP = zeros(Int64,O,2)
    ΔG = zeros(O)
    # Preallocate vector to store chemical potential at standard conditions
    μ0 = zeros(M)
    # Set fixed intervals
    rs = 0.5*ones(M-1)
    # Take maximum value to be μrange, this is unphysical but works as we are intrested in changes
    μ0[1] = μrange
    for i = 2:M
        μ0[i] = μ0[i-1] - rs[i-1]*μrange
    end
    # Now deterministically assign reactions
    for i = 1:2
        RP[i,1] = i
        RP[i,2] = i+1
        ΔG[i] = μ0[i+1] - μ0[i]
    end
    return(RP,ΔG)
end

# similar to the other η but now with a substrate consuming species
function initialise_η2(N::Int64,mq::Float64,mK::Float64,mk::Float64)
    # Three metabolities alterable by two reactions in this case
    M = 3
    O = 2
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # And all proportionality constants are the same for simplicity
    g = ones(N)
    # All but resource 1 is not supplied
    κ = zeros(M)
    κ[1] = 100.0
    # Assume that all δ's are equal
    δi = 1.0
    δ = δi*ones(M)
    # Use deterministic maintenance rates
    mm = 1.0
    m = mm*ones(N)
    # Generate random set of reactions
    μrange = 6e5 # Slightly larger μrange so that two cases can be compared
    RP, ΔG = two_reactions(μrange)
    # Preallocate vector of reactions
    reacs = Array{Reaction,1}(undef,O)
    for i = 1:O
        reacs[i] = make_Reaction(i,RP[i,1],RP[i,2],ΔG[i])
    end
    # Preallocate vector of microbes
    mics = Array{Microbe,1}(undef,N)
    # Then construct microbes
    for i = 1:N
        # Only one reaction which each microbe uses
        R = 1
        if i == 1
            Reacs = [2]
        else
            Reacs = [1]
        end
        # Kinetic parameters are fixed
        qm = [mq]
        KS = [mK]
        kr = [mk]
        # Find corresponding η's for these reactions
        if i == 1
            η = [3.5]
        else
            η = [single_ηs(reacs,T,N-1,i-1)]
        end
        # Can finally generate microbe
        mics[i] = make_Microbe(m[i],g[i],R,Reacs,η,qm,KS,kr)
    end
    # Now make the parameter set
    ps = make_InhibParameters(N,M,O,T,κ,δ,reacs,mics)
    return(ps)
end
