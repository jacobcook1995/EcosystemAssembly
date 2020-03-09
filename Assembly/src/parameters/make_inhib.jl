# Script that makes the parameters needed for different types of simulation
# of the model with inhibition
using Assembly

export initialise

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

# function to run simulation of the Marsland model
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
