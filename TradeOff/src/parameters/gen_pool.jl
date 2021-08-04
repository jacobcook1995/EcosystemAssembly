# Script to generate and save a pool of microbes with particular parameter ranges

export new_pool, new_mic, fix_reactions

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

# function that chooses uniformly mixed values for η
function choose_η_mix(reacs::Array{Reaction,1},Reacs::Array{Int64,1},T::Float64,mratio::Float64)
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
        # And then use to determine an upper bound on η
        ηh = -(dG + Rgas*T*log(mratio))/(ΔGATP)
        η[i] = (ηh-ηl)*rand() + ηl
    end
    return(η)
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

# Function to generate a new pool
function new_pool(Nt::Int64,M::Int64,Rs::Array{Int64,1},d::Float64,μrange::Float64,mratio::Float64)
    # First generate random unique indetifier for this pool
    PID = randstring(['0':'9'; 'a':'f'])
    # Print out that this is happening
    println("Generating random pool with identifer: $(PID)")
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
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Number of doublings required to dilute to 1%
    fd = log(100)/log(2)
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # The number of ATP per translation step, including the cost of amino acid sythesis
    # This figure is taken from Lynch and Marinov 2015
    χl = 29.0
    # Estimate maximum additional cost to be three times mimumum (from Fig 1 in Roller et al)
    χu = 3.5*χl
    # This is a slightly arbitary choice for Kγ
    Kγ = 5e8
    # Set mimumum KΩ value
    KΩm = 1e9
    # Use formula to calculate how many reactions are implied
    O = 2*M - 3
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # Generate fixed set of reactions
    RP, ΔG = fix_reactions(O,M,μrange,T)
    # Preallocate vector of reactions
    reacs = Array{Reaction,1}(undef,O)
    for i = 1:O
        reacs[i] = make_Reaction(i,RP[i,1],RP[i,2],ΔG[i])
    end
    # Assume that half saturation occurs at a quarter κ/δ
    KS = (1/4)*5.5e-3
    # Arbitary number that seems to give decent survival
    kc = 10.0
    # The reversibility factor remains the same as previously
    kr = 10.0
    # Preallocate vector of microbes
    mics = Array{Microbe,1}(undef,Nt)
    # Then construct microbes
    for i = 1:Nt
        # For each microbe generate random set of reactions
        R, Reacs = choose_reactions(O,Rs)
        # Roller et al suggest a factor of 10 variation in max growth rate
        # So varies between 0.1 and 1.0
        ω = 0.1 + 0.9*rand()
        # Want KΩ variations to be smaller, i.e. factor of 2.5
        KΩ = 1.5*(ω-0.1)*KΩm + 1*KΩm
        # Make vectors of the (fixed) kinetic parameters
        kcs, KSs, krs = kin_rand(kc,KS,kr,R)
        # Reactions given random proportional weightings, done this in the simplest way possible
        ϕP = rand(R)
        ϕP = ϕP/sum(ϕP)
        # Find corresponding η's for these reactions
        η = choose_η_mix(reacs,Reacs,T,mratio)
        # Can finally generate microbe
        mics[i] = make_Microbe(MC,γm,Kγ,χl,χu,Pb,d,ϕH,KΩ,fd,ω,R,Reacs,η,kcs,KSs,krs,n,ϕP,i,PID)
    end
    # Write out necessary data
    jldopen("Pools/ID=$(PID)N=$(Nt)M=$(M)d=$(d)u=$(μrange).jld","w") do file
        # Write out species pool
        write(file,"mics",mics)
        # And the metabolic network details for later verification
        write(file,"M",M)
        write(file,"reacs",reacs)
    end
    return(nothing)
end

# Alternative function to generate a single microbe
function new_mic(M::Int64,Rs::Array{Int64,1},d::Float64,μrange::Float64,mratio::Float64)
    # First generate random unique indetifier for this pool
    PID = randstring(['0':'9'; 'a':'f'])
    # Print out that this is happening
    println("Generating random pool with identifer: $(PID)")
    # Only generating one microbe so will be ID: 1
    ID = 1
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
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Number of doublings required to dilute to 1%
    fd = log(100)/log(2)
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # The number of ATP per translation step, including the cost of amino acid sythesis
    # This figure is taken from Lynch and Marinov 2015
    χl = 29.0
    # Estimate maximum additional cost to be three times mimumum (from Fig 1 in Roller et al)
    χu = 3.5*χl
    # This is a slightly arbitary choice for Kγ
    Kγ = 5e8
    # Set mimumum KΩ value
    KΩm = 1e9
    # Use formula to calculate how many reactions are implied
    O = 2*M - 3
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # Generate fixed set of reactions
    RP, ΔG = fix_reactions(O,M,μrange,T)
    # Preallocate vector of reactions
    reacs = Array{Reaction,1}(undef,O)
    for i = 1:O
        reacs[i] = make_Reaction(i,RP[i,1],RP[i,2],ΔG[i])
    end
    # Assume that half saturation occurs at a quarter κ/δ
    KS = (1/4)*5.5e-3
    # Arbitary number that seems to give decent survival
    kc = 10.0
    # The reversibility factor remains the same as previously
    kr = 10.0
    # Generate random set of reactions
    R, Reacs = choose_reactions(O,Rs)
    # Roller et al suggest a factor of 10 variation in max growth rate
    # So varies between 0.1 and 1.0
    ω = 0.1 + 0.9*rand()
    # Want KΩ variations to be smaller, i.e. factor of 2.5
    KΩ = 1.5*(ω-0.1)*KΩm + 1*KΩm
    # Make vectors of the (fixed) kinetic parameters
    kcs, KSs, krs = kin_rand(kc,KS,kr,R)
    # Reactions given random proportional weightings, done this in the simplest way possible
    ϕP = rand(R)
    ϕP = ϕP/sum(ϕP)
    # Find corresponding η's for these reactions
    η = choose_η_mix(reacs,Reacs,T,mratio)
    # Can finally generate microbe
    mic = make_Microbe(MC,γm,Kγ,χl,χu,Pb,d,ϕH,KΩ,fd,ω,R,Reacs,η,kcs,KSs,krs,n,ϕP,ID,PID)
    return(mic)
end
