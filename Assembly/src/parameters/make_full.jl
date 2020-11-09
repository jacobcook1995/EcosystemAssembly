# Script that makes the parameters needed for simulation of the full proteome model
export initialise

# function to randomly choose the reactions that a specific microbe can make use of
 function choose_reactions(O::Int64,Rl::Int64,Ru::Int64)
     @assert Ru <= O "Strain cannot have more reactions than the number of possible reactions"
     # Make required Uniform distribution between limits Rl and Ru
     R = rand(Rl:Ru)
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

# function that chooses uniformly mixed values for η
function choose_η_mix(reacs::Array{Reaction,1},Reacs::Array{Int64,1},T::Float64,syn::Bool)
    # Preallocate memory to store η's
    η = zeros(length(Reacs))
    # Set a constant lower bound
    ηl = 1/3
    # Placeholder value
    mratio = 0
    # Choose actual value based on whether we want syntrophy on or not
    if syn == true
        # Set minimum equilibrium product to substrate ratio
        mratio = 1e-2
    else
        mratio = 1e5
    end
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

# function to take in average kinetic parameters and return randomised vectors of them
function kin_rand(kc::Float64,KS::Float64,kr::Float64,R::Int64)
    # Make probability distribution
    d = Normal()
    # Generate array of random numbers
    rs = rand(d,R,3)
    # Then use to generate randomly varibles that have same probability to be half as double average
    kcs = kc*(2.0.^rs[:,1])
    KSs = KS*(2.0.^rs[:,2])
    krs = kr*(2.0.^rs[:,3])
    return(kcs,KSs,krs)
end

# function to generate parameter set for the model with inhibition
function initialise(N::Int64,M::Int64,O::Int64,Rl::Int64,Ru::Int64,kc::Float64,KS::Float64,kr::Float64,syn::Bool)
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
    # Taken from my ATP data, we'll see how this goes
    Kγ = 4.90e7
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # Taken from my fits to ATP data, we'll see how this goes
    KΩ = 1.83e7
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
        R, Reacs = choose_reactions(O,Rl,Ru)
        # Make vectors of the (fixed) kinetic parameters
        kcs, KSs, krs = kin_rand(kc,KS,kr,R)
        # Reactions given random proportional weightings, done this in the simplest way possible
        ϕP = rand(R)
        ϕP = ϕP/sum(ϕP)
        # Find corresponding η's for these reactions
        η = choose_η_mix(reacs,Reacs,T,syn)
        # Can finally generate microbe
        mics[i] = make_MicrobeP(MC,γm,ρ,Kγ,Pb,d,ϕH,KΩ,fd,R,Reacs,η,kcs,KSs,krs,n,ϕP)
    end
    # Now make the parameter set
    ps = make_FullParameters(N,M,O,T,κ,δ,reacs,mics)
    return(ps)
end
