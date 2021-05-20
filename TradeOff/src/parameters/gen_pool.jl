# Script to generate and save a pool of microbes with particular parameter ranges

export new_pool

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

# HUGE AMOUNT OF STUFF HARDCODED IN HERE AT THE MOMENT, THIS WILL HAVE TO CHANGE
# Function to generate a new pool
function new_pool(Nt::Int64,M::Int64,Rl::Int64,Ru::Int64)
    # First generate random unique indetifier for this pool
    PID = randstring(['0':'9'; 'a':'f'])
    # Print out that this is happening
    println("Generating random pool with identifer: $(PID)")
    # Assume that half saturation occurs at a quarter κ/δ
    KS = (1/4)*5.5e-3
    # Arbitary number that seems to give decent survival
    kc = 10.0
    # The reversibility factor remains the same as previously
    kr = 10.0
    # Always want to allow near to equilbrium reactions now
    syn = true
    # Use formula to calculate how many reactions are implied
    O = 2*M - 3
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
    # Back to old arbitary figures, think this might be the best route
    Kγ = 5e8
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # Back to old arbitary figures, think this might be the best route
    KΩ = 1e9
    # Number of doublings required to dilute to 1%
    fd = log(100)/log(2)
    # Chosen so that 100 steps yields slightly more free energy than respiring glucose
    μrange = 5e6*(M/25)
    # Generate fixed set of reactions
    RP, ΔG = fix_reactions(O,M,μrange,T)
    # Preallocate vector of reactions
    reacs = Array{Reaction,1}(undef,O)
    for i = 1:O
        reacs[i] = make_Reaction(i,RP[i,1],RP[i,2],ΔG[i])
    end
    # Preallocate vector of microbes
    mics = Array{Microbe,1}(undef,Nt)
    # Then construct microbes
    for i = 1:Nt
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
        mics[i] = make_Microbe(MC,γm,ρ,Kγ,Pb,d,ϕH,KΩ,fd,R,Reacs,η,kcs,KSs,krs,n,ϕP,i,PID)
    end
    # Write out necessary data
    jldopen("Pools/N=$(Nt)M=$(M)Reacs$(Rl)-$(Ru)ID=$(PID).jld","w") do file
        # Write out species pool
        write(file,"mics",mics)
        # And the metabolic network details for later verification
        write(file,"M",M)
        write(file,"reacs",reacs)
    end
    return(nothing)
end
