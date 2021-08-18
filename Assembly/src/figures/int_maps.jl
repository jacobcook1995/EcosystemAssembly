# Script to make the interaction heat maps
using Assembly
using JLD
using Plots
using Plots.PlotMeasures
import PyPlot

# function to make heat maps of interactions between functonal groups
function interaction_maps()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 6
        error("Insufficent inputs provided (looking for 6)")
    end
    # Preallocate the variables I want to extract from the input
    Rl = 0
    Ru = 0
    syn = true
    rps = 0
    Ni = 0
    en = ARGS[6]
    # Check that all arguments can be converted to integers
    try
        Rl = parse(Int64,ARGS[1])
        Ru = parse(Int64,ARGS[2])
        syn = parse(Bool,ARGS[3])
        rps = parse(Int64,ARGS[4])
        Ni = parse(Int64,ARGS[5])
    catch e
           error("need to provide 4 integers, a bool and a string")
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
    # Check that initial number of strains is sensible
    if Ni < 1
        error("number of strains cannot be less than 1")
    end
    println("Compiled!")
    flush(stdout)
    # Load first parameter file to find relevant info
    pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run1Ns$(Ni).jld"
    if ~isfile(pfile)
        error("run 1 is missing a parameter file")
    end
    ps = load(pfile,"ps")
    # Preallocate data storage
    comps = zeros(ps.M-1,ps.M-1)
    facls = zeros(ps.M-1,ps.M-1)
    synps = zeros(ps.M-1,ps.M-1)
    polls = zeros(ps.M-1,ps.M-1)
    no_comps = zeros(ps.M-1,ps.M-1)
    no_facls = zeros(ps.M-1,ps.M-1)
    no_synps = zeros(ps.M-1,ps.M-1)
    no_polls = zeros(ps.M-1,ps.M-1)
    # Loop over the repeats to extract neccessary data
    for i = 1:rps
        # Find relevant file of interactions
        ifile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/IntsReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
        if ~isfile(pfile)
            error("run $(i) is missing an interactions file")
        end
        # Read in relevant files
        pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        # Load in interactions and interaction strengths
        ints = load(ifile,"ints")
        in_str = load(ifile,"in_str")
        ps = load(pfile,"ps")
        # Loop over every strain
        for j = 1:ps.N
            # Find the ith strain's main reaction
            _, ind = findmax(ps.mics[j].ϕP)
            # And use this to find it's functional group
            fg_j = ps.reacs[ps.mics[j].Reacs[ind]].Rct
            # Loop over all other strains
            for k = 1:ps.N
                # Find the jth strain's functional group
                _, ind = findmax(ps.mics[k].ϕP)
                # And use this to find it's functional group
                fg_k = ps.reacs[ps.mics[k].Reacs[ind]].Rct
                # Loop over all metabolites that the function group can interact via
                for l = 1:ps.M
                    # Check interaction types (no interaction indicated by 0)
                    if ints[j,k,l] == 1
                        # Add interaction strength to competition container
                        comps[fg_j,fg_k] += in_str[j,k,l]
                        no_comps[fg_j,fg_k] += 1
                    elseif ints[j,k,l] == 2
                        # Add interaction strength to facilitation container
                        facls[fg_j,fg_k] += in_str[j,k,l]
                        no_facls[fg_j,fg_k] += 1
                    elseif ints[j,k,l] == 3
                        # Add interaction strength to syntrophy container
                        synps[fg_j,fg_k] += in_str[j,k,l]
                        no_synps[fg_j,fg_k] += 1
                    elseif ints[j,k,l] == 4
                        # Add interaction strength to pollution container
                        polls[fg_j,fg_k] += in_str[j,k,l]
                        no_polls[fg_j,fg_k] += 1
                    end
                end
            end
        end
    end
    # Now setup plotting
    pyplot(dpi=200)
    # Set up the range to use for the axes (only showing first 12)
    rge = 1:12
    # Make competition plot
    p1 = heatmap(rge,rge,comps[12:-1:1,1:12],c=:heat,colorbar_title="Total interaction strength",title="Competiton")
    plot!(p1,xticks=(2:2:12),yticks=(2:2:12,12:-2:2),xlabel="Functional group",ylabel="Functional group impacted")
    # Add annotation
    annotate!(p1,-0.5,13,text("A",17,:black))
    # Make facilitation plot
    p2 = heatmap(rge,rge,facls[12:-1:1,1:12],c=:heat,colorbar_title="Total interaction strength",title="Facilitation")
    plot!(p2,xticks=(2:2:12),yticks=(2:2:12,12:-2:2),xlabel="Functional group",ylabel="Functional group impacted")
    # Add annotation
    annotate!(p2,-0.5,13,text("B",17,:black))
    # Make syntrophy plot
    p3 = heatmap(rge,rge,synps[12:-1:1,1:12],c=:heat,colorbar_title="Total interaction strength",title="Syntrophy")
    plot!(p3,xticks=(2:2:12),yticks=(2:2:12,12:-2:2),xlabel="Functional group",ylabel="Functional group impacted")
    # Add annotation
    annotate!(p3,-0.5,13,text("C",17,:black))
    # Make pollution plot
    p4 = heatmap(rge,rge,polls[12:-1:1,1:12],c=:heat,colorbar_title="Total interaction strength",title="Pollution")
    plot!(p4,xticks=(2:2:12),yticks=(2:2:12,12:-2:2),xlabel="Functional group",ylabel="Functional group impacted")
    # Add annotation
    annotate!(p4,-0.5,13,text("D",17,:black))
    # Combine all figures into one
    pc1 = plot(p1,p2,p3,p4,layout=4,size=(1200,800),margin=7.5mm)
    savefig(pc1,"Output/SI/heat_strength.png")
    # Do the same load of heatmaps but for the number of interactions, starting with competition
    p5 = heatmap(rge,rge,no_comps[12:-1:1,1:12],c=:heat,colorbar_title="Total number of interactions",title="Competiton")
    plot!(p5,xticks=(2:2:12),yticks=(2:2:12,12:-2:2),xlabel="Functional group",ylabel="Functional group impacted")
    # Add annotation
    annotate!(p5,-0.5,13,text("A",17,:black))
    # Make facilitation plot
    p6 = heatmap(rge,rge,no_facls[12:-1:1,1:12],c=:heat,colorbar_title="Total number of interactions",title="Facilitation")
    plot!(p6,xticks=(2:2:12),yticks=(2:2:12,12:-2:2),xlabel="Functional group",ylabel="Functional group impacted")
    # Add annotation
    annotate!(p6,-0.5,13,text("B",17,:black))
    # Make syntrophy plot
    p7 = heatmap(rge,rge,no_synps[12:-1:1,1:12],c=:heat,colorbar_title="Total number of interactions",title="Syntrophy")
    plot!(p7,xticks=(2:2:12),yticks=(2:2:12,12:-2:2),xlabel="Functional group",ylabel="Functional group impacted")
    # Add annotation
    annotate!(p7,-0.5,13,text("C",17,:black))
    # Make pollution plot
    p8 = heatmap(rge,rge,no_polls[12:-1:1,1:12],c=:heat,colorbar_title="Total number of interactions",title="Pollution")
    plot!(p8,xticks=(2:2:12),yticks=(2:2:12,12:-2:2),xlabel="Functional group",ylabel="Functional group impacted")
    # Add annotation
    annotate!(p8,-0.5,13,text("D",17,:black))
    # Combine all figures into one
    pc2 = plot(p5,p6,p7,p8,layout=4,size=(1200,800),margin=7.5mm)
    savefig(pc2,"Output/SI/heat_number.png")
    return(nothing)
end

@time interaction_maps()
