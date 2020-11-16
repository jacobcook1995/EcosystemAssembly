# Script to find and analyse the interaction strength and type of strains in my simulation
using Assembly
using JLD
using SymPy
using LinearAlgebra

# function to quantify the interactions
function quantify_ints()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 4
        error("Insufficent inputs provided (looking for 4)")
    end
    # Preallocate the variables I want to extract from the input
    Rl = 0
    Ru = 0
    syn = true
    rps = 0
    # Check that all arguments can be converted to integers
    try
        Rl = parse(Int64,ARGS[1])
        Ru = parse(Int64,ARGS[2])
        syn = parse(Bool,ARGS[3])
        rps = parse(Int64,ARGS[4])
    catch e
           error("need to provide 3 integers and a bool")
    end
    # Check that simulation type is valid
    if Rl < 1
        error("lower bound has to be at least one reaction")
    end
    if Ru < Rl
        error("upper bound can't be lower than the lower bound")
    end
    # Check that number of simulations is greater than 0
    if rps < 1
        error("number of repeats can't be less than 1")
    end
    println("Compiled!")
    # Loop over the repeats
    for i = 1:rps
        # Read in relevant files
        pfile = "Data/$(Rl)-$(Ru)$(syn)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        ofile = "Data/$(Rl)-$(Ru)$(syn)/RedOutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        if ~isfile(ofile)
            error("run $(i) is missing an output file")
        end
        efile = "Data/$(Rl)-$(Ru)$(syn)/RedExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        if ~isfile(efile)
            error("run $(i) is missing an extinct file")
        end
        # Basically just loading everything out as I'm not sure what I'll need
        ps = load(pfile,"ps")
        C = load(ofile,"C")
        T = load(ofile,"T")
        out = load(ofile,"out")
        inf_out = load(ofile,"inf_out")
        ded = load(efile,"ded")
        # Preallocate vector of forces
        F = Array{Sym,1}(undef,3*ps.N+ps.M)
        # Now find vector of forces
        F = Force(ps,F)
        # Preallocate perturbation
        pout = zeros(3*ps.N+ps.M)
        # Populations, energies, and fractions remain fixed
        pout[1:ps.N] = inf_out[1:ps.N]
        pout[(ps.N+ps.M+1):(2*ps.N+ps.M)] = inf_out[(ps.N+ps.M+1):(2*ps.N+ps.M)]
        pout[(2*ps.N+ps.M+1):end] = inf_out[(2*ps.N+ps.M+1):end]
        # Preallocate array to store interaction types
        ints = zeros(Int64,ps.N,ps.N,ps.M)
        # Preallocate array to store interaction strengths
        in_str = zeros(ps.N,ps.N,ps.M)
        # Set fraction to increase concentrations by
        frc = 0.1
        # Loop over metabolites
        for j = 1:ps.M
            # Skip metabolites with zero concentration
            if inf_out[ps.N+j] != 0.0
                # Loop over metabolites to make perturbation
                for k = 1:ps.M
                    if k != j
                        pout[ps.N+k] = inf_out[ps.N+k]
                    else
                        pout[ps.N+k] = inf_out[ps.N+k]*(1.0+frc)
                    end
                end
                # Find perturbed forces
                nF = nForce(F,pout,ps)
                # Only interested in the forces on the ATP
                Fatp = nF[(ps.N+ps.M+1):(2*ps.N+ps.M)]
                # Preallocate postive/negative/zero vector
                pn0 = zeros(Int64,ps.N)
                # Preallocate supply/consumption/zero vector
                sc0 = zeros(Int64,ps.N)
                # Preallocate fluxes
                fl = zeros(ps.N)
                # Loop over the strains
                for k = 1:ps.N
                    # Check whether the effect on the strain is +ve, -ve or ~0
                    if Fatp[k] > 1e-5
                        # Positive denoted by 1
                        pn0[k] = 1
                    elseif Fatp[k] < -1e-5
                        # Negative denoted by -1
                        pn0[k] = -1
                    end
                    # Check if strain interacts with metabolite
                    if any((ps.reacs[ps.mics[k].Reacs].↦:Rct) .== j) || any((ps.reacs[ps.mics[k].Reacs].↦:Prd) .== j)
                        # Sum rate terms to find net effect
                        qT = 0
                        # If it does then loop over reactions
                        for l = 1:ps.mics[k].R
                            # Find amount of enzyme
                            E = Eα(inf_out[2*ps.N+ps.M+k],ps.mics[k],l)
                            # Find reaction
                            r = ps.reacs[ps.mics[k].Reacs[l]]
                            # Find reaction rate
                            q = qs(inf_out[ps.N+r.Rct],inf_out[ps.N+r.Prd],E,l,ps.mics[k],ps.T,r)
                            # Check if reactant or product
                            if r.Rct == j
                                qT -= q
                            elseif r.Prd == j
                                qT += q
                            end
                        end
                        # Check if there's net supply or consumption
                        fl[k] = qT
                        if qT > 0.0
                            # Supply indicated by 1
                            sc0[k] = 1
                        elseif qT < 0.0
                            # And consumption by -1
                            sc0[k] = -1
                        end
                    end
                end
                # Loop over vectors to make interaction matrix
                for k = 1:ps.N
                    for l = 1:ps.N
                        # Save the interaction strength
                        in_str[k,l,j] = abs(Fatp[k])*(abs(fl[l])/sum(abs.(fl)))
                        # Check for positive cases
                        if pn0[k] == 1
                            # Crossfeeding
                            if sc0[l] == 1
                                ints[k,l,j] = 2
                            # Competition
                            elseif sc0[l] == -1
                                ints[k,l,j] = 1
                            end
                        # And negative case
                        elseif pn0[k] == -1
                            # Competiton via pollution
                            if sc0[l] == 1
                                ints[k,l,j] = 4
                            # Obligate syntrophy
                            elseif sc0[l] == -1
                                ints[k,l,j] = 3
                            end
                        end
                    end
                end
            end
        end
        println(ints)
    end
    return(nothing)
end


@time quantify_ints()
