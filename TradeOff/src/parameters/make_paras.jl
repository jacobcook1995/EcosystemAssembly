# Script that makes the parameters needed for simulation of the full proteome model
export initialise

# THIS PROCEDURE WILL HAVE TO CHANGE WHEN I INCLUDE NEW PARAMETERS

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

# function to generate parameter set for the fixed parameters
function initialise(M::Int64,O::Int64)
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # From Posfai et al (2017) dilution rate 0.21 per hour
    δ = 6.0e-5*ones(M) # Metabolite dilution rate
    # Human blood glucose is approximatly 5.5 m mol per litre (wikipedia)
    # Sensible order of magnitude to aim for set κ/δ = 5.5e-3
    κ = zeros(M)
    # All but resource 1 is not supplied
    κ[1] = 3.3e-7 # Metabolite supply rate
    # Chosen so that 100 steps yields slightly more free energy than respiring glucose
    μrange = 5e6*(M/25)
    # Generate fixed set of reactions
    RP, ΔG = fix_reactions(O,M,μrange,T)
    # Preallocate vector of reactions
    reacs = Array{Reaction,1}(undef,O)
    for i = 1:O
        reacs[i] = make_Reaction(i,RP[i,1],RP[i,2],ΔG[i])
    end
    # Now make the parameter set
    ps = make_TOParameters(M,O,T,κ,δ,reacs)
    return(ps)
end

# function to generate parameter set for the fixed parameters
function initialise(M::Int64,O::Int64,μrange::Float64)
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # From Posfai et al (2017) dilution rate 0.21 per hour
    δ = 6.0e-5*ones(M) # Metabolite dilution rate
    # Human blood glucose is approximatly 5.5 m mol per litre (wikipedia)
    # Sensible order of magnitude to aim for set κ/δ = 5.5e-3
    κ = zeros(M)
    # All but resource 1 is not supplied
    κ[1] = 3.3e-7 # Metabolite supply rate
    # Generate fixed set of reactions
    RP, ΔG = fix_reactions(O,M,μrange,T)
    # Preallocate vector of reactions
    reacs = Array{Reaction,1}(undef,O)
    for i = 1:O
        reacs[i] = make_Reaction(i,RP[i,1],RP[i,2],ΔG[i])
    end
    # Now make the parameter set
    ps = make_TOParameters(M,O,T,κ,δ,reacs)
    return(ps)
end
