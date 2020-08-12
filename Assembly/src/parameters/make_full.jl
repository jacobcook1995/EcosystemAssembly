# Script that makes the parameters needed for simulation of the full proteome model
export initialise

# Function that chooses uniformly low values for η
function choose_η_low(reacs::Array{Reaction,1},Reacs::Array{Int64,1},T::Float64)
    # Preallocate memory to store η's
    η = zeros(length(Reacs))
    # Set a constant lower bound
    ηl = 1/3
    # Loop over all reactions
    for i = 1:length(η)
        # Identify which reaction we are considering
        I = Reacs[i]
        # Find corresponding Gibbs free energy change
        dG = reacs[I].ΔG0
        ηh = 0.4*(-dG/ΔGATP)
        # Considering uniformly distributed η's'
        η[i] = (ηh-ηl)*rand() + ηl
    end
    return(η)
end

# function that chooses uniformly moderate values for η
function choose_η_med(reacs::Array{Reaction,1},Reacs::Array{Int64,1},T::Float64)
    # Preallocate memory to store η's
    η = zeros(length(Reacs))
    # Loop over all reactions
    for i = 1:length(η)
        # Identify which reaction we are considering
        I = Reacs[i]
        # Find corresponding Gibbs free energy change
        dG = reacs[I].ΔG0
        # Use to find upper and lower bounds
        ηl = 0.4*(-dG/ΔGATP)
        ηh = 0.7*(-dG/ΔGATP)
        # Considering uniformly distributed η's'
        η[i] = (ηh-ηl)*rand() + ηl
    end
    return(η)
end

# function that chooses uniformly high values for η
function choose_η_high(reacs::Array{Reaction,1},Reacs::Array{Int64,1},T::Float64)
    # Preallocate memory to store η's
    η = zeros(length(Reacs))
    # Set minimum equilibrium product to substrate ratio
    mratio = 1e-2
    # Loop over all reactions
    for i = 1:length(η)
        # Identify which reaction we are considering
        I = Reacs[i]
        # Find corresponding Gibbs free energy change
        dG = reacs[I].ΔG0
        # Use to find lower bound
        ηl = 0.7*(-dG/ΔGATP)
        # And then use to determine an upper bound on η
        ηh = -(dG + Rgas*T*log(mratio))/(ΔGATP)
        η[i] = (ηh-ηl)*rand() + ηl
    end
    return(η)
end

# function that chooses uniformly mixed values for η
function choose_η_mix(reacs::Array{Reaction,1},Reacs::Array{Int64,1},T::Float64)
    # Preallocate memory to store η's
    η = zeros(length(Reacs))
    # Set a constant lower bound
    ηl = 1/3
    # Set minimum equilibrium product to substrate ratio
    mratio = 1e-2
    # Loop over all reactions
    for i = 1:length(η)
        # Identify which reaction we are considering
        I = Reacs[i]
        # Find corresponding Gibbs free energy change
        dG = reacs[I].ΔG0
        # And then use to determine an upper bound on η
        ηh = -(dG + Rgas*T*log(mratio))/(ΔGATP)
        η[i] = (ηh-ηl)*rand() + ηl
    end
    return(η)
end

# This one is just a toy function to generate fixed low values of η
function choose_η_fix(reacs::Array{Reaction,1},Reacs::Array{Int64,1},T::Float64)
    # Preallocate memory to store η's
    η = zeros(length(Reacs))
    for i = 1:length(η)
        # Identify which reaction we are considering
        I = Reacs[i]
        # Find corresponding Gibbs free energy change
        dG = reacs[I].ΔG0
        # And use to determine the value of η
        η[i] = 0.9*(-dG/ΔGATP)
    end
    return(η)
end

# function to generate fix set of reaction for our model. Each metabolite can be broken
# down into metabolites 1 or 2 steps down. The steps between metabolites are fixed.
function fix_reactions(O::Int64,M::Int64,μrange::Float64,T::Float64)
    @assert O == 2*M - 3 "Miscalulated the number of reactions expected"
    # preallocate output
    RP = zeros(Int64,O,2)
    ΔG = zeros(O)
    # find ΔG for a single step
    dG = -μrange/(M-1)
    # loop over all reactions
    for i = 1:O
        # find odd numbers
        if i % 2 != 0
            RP[i,1] = ceil(i/2)
            RP[i,2] = ceil(i/2) + 1
            ΔG[i] = dG
        else
            RP[i,1] = ceil(i/2)
            RP[i,2] = ceil(i/2) + 2
            ΔG[i] = 2*dG
        end
    end
    return(RP,ΔG)
end

# function to generate parameter set for the model with inhibition
function initialise(N::Int64,M::Int64,O::Int64,mR::Float64,sdR::Float64,kc::Float64,KS::Float64,kr::Float64,choice::Int64)
    @assert O >= mR + 5*sdR "Not enough reactions to ensure that microbes have on average mR reactions"
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # Cell mass is taken from Bremer H, Dennis P (1996) Modulation of chemical
    # composition and other parameters of the cell by growth rate (Book chapter).
    MC = 10^8
    # Elongation rate also taken from Bremer (1996), convert from minutes to seconds
    γm = 1260.0/60.0
    # Now preallocate protein masses
    n = zeros(Int64,3)
    # ribosome mass taken from Keseler IM, et al. (2011)
    n[1] = 7459
    # Other protein mass averaged from Brandt F, et al. (2009)
    n[2:3] .= 300
    # Need to think carefully about what a reasonable parameter value is here
    # This has a large impact on the fraction, so should be careful with this one
    d = 6.0e-5
    # The number of ATP per translation step, including the cost of amino acid sythesis
    # This figure is taken from Lynch and Marinov 2015
    ρ = 29.0
    # This is currently a paramter which I am fiddling
    Kγ = 5e8
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # Give omega a fairly arbitary value for now, would be expected to be of similar order to Kγ
    KΩ = 2*Kγ
    # Number of doublings required to dilute to 1%
    fd = log(100)/log(2)
    # From Posfai et al (2017) dilution rate 0.21 per hour
    δ = 6.0e-5*ones(M) # Metabolite dilution rate
    # Human blood glucose is approximatly 5.5 m mol per litre (wikipedia)
    # Sensible order of magnitude to aim for set κ/δ = 5.5e-3
    κ = zeros(M)
    # All but resource 1 is not supplied
    κ[1] = 3.3e-7 # Metabolite supply rate
    # Chosen so that 100 steps yields slightly more free energy than respiring glucose
    μrange = 3e6*(M/25)
    # Generate fixed set of reactions
    RP, ΔG = fix_reactions(O,M,μrange,T)
    # Preallocate vector of reactions
    reacs = Array{Reaction,1}(undef,O)
    for i = 1:O
        reacs[i] = make_Reaction(i,RP[i,1],RP[i,2],ΔG[i])
    end
    # Preallocate vector of microbes
    mics = Array{MicrobeP,1}(undef,N)
    # Then construct microbes
    for i = 1:N
        # For each microbe generate random set of reactions
        R, Reacs = choose_reactions(O,mR,sdR)
        # Make vectors of the (fixed) kinetic parameters
        kcs = kc*ones(R)
        KSs = KS*ones(R)
        krs = kr*ones(R)
        # Reactions given random proportional weightings, done this in the simplest way possible
        ϕP = rand(R)
        ϕP = ϕP/sum(ϕP)
        # Find corresponding η's for these reactions
        if choice == 1
            η = choose_η_low(reacs,Reacs,T)
        elseif choice == 2
            η = choose_η_med(reacs,Reacs,T)
        elseif choice == 3
            η = choose_η_high(reacs,Reacs,T)
        else
            η = choose_η_mix(reacs,Reacs,T)
        end
        # Can finally generate microbe
        mics[i] = make_MicrobeP(MC,γm,ρ,Kγ,Pb,d,ϕH,KΩ,fd,R,Reacs,η,kcs,KSs,krs,n,ϕP)
    end
    # Now make the parameter set
    ps = make_FullParameters(N,M,O,T,κ,δ,reacs,mics)
    return(ps)
end