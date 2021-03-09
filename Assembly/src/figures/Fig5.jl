# Script to construct figure 5
using Assembly
using Plots
using JLD
using LaTeXStrings
using Statistics
import PyPlot

function figure5(Rl::Int64,Ru::Int64,syn::Bool,rps::Int64,Ni::Int64,en::String)
    println("Compiled!")
    # Preallocate memory to store number of interactions
    ins1 = zeros(rps)
    ins2 = zeros(rps)
    ins3 = zeros(rps)
    ins4 = zeros(rps)
    # Preallocate memory to store mean interaction strengths
    mean1 = fill(NaN,rps)
    mean2 = fill(NaN,rps)
    mean3 = fill(NaN,rps)
    mean4 = fill(NaN,rps)
    # Preallocate vectors to store all interaction strengths
    sts1 = []
    sts2 = []
    sts3 = []
    sts4 = []
    # Loop over parameter sets
    for i = 1:rps
        # Read in relevant files
        ifile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/IntsReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
        if ~isfile(ifile)
            error("run $(i) is missing an interaction file")
        end
        # Just need to load ints and their strengths
        ints = load(ifile,"ints")
        in_str = load(ifile,"in_str")
        # Store number of each interaction type
        ins1[i] = count(x->x==1,ints)
        ins2[i] = count(x->x==2,ints)
        ins3[i] = count(x->x==3,ints)
        ins4[i] = count(x->x==4,ints)
        # Find corresponding positions
        pos1 = findall(x->x==1,ints)
        pos2 = findall(x->x==2,ints)
        pos3 = findall(x->x==3,ints)
        pos4 = findall(x->x==4,ints)
        # Use to find vector of interaction strengths
        str1 = in_str[pos1]
        str2 = in_str[pos2]
        str3 = in_str[pos3]
        str4 = in_str[pos4]
        # Check if there are any interactions of that type
        if length(str1) >= 1
            # Find mean strength of each interaction
            mean1[i] = mean(str1)
            sts1 = cat(sts1,str1,dims=1)
        end
        # Do for each interaction type
        if length(str2) >= 1
            mean2[i] = mean(str2)
            sts2 = cat(sts2,str2,dims=1)
        end
        if length(str3) >= 1
            mean3[i] = mean(str3)
            sts3 = cat(sts3,str3,dims=1)
        end
        if length(str4) >= 1
            mean4[i] = mean(str4)
            sts4 = cat(sts4,str4,dims=1)
        end
    end
    # Setup plotting
    pyplot()
    theme(:wong2,dpi=200)
    # Make plot title
    tl = ""
    if syn == true
        tl = "$(Rl)-$(Ru) reactions per strain"
    else
        tl = "$(Rl)-$(Ru) reactions per strain (no syntrophy)"
    end
    # Now plot both interactions types and strengths
    plot(title=tl,xlabel="Number of interactions",ylabel="Number of ecosystems")
    histogram!(ins1,fillalpha=0.75,label="Competiton")
    histogram!(ins2,fillalpha=0.75,label="Facilitation")
    histogram!(ins3,fillalpha=0.75,label="Syntrophy")
    histogram!(ins4,fillalpha=0.75,label="Pollution")
    savefig("Output/Fig5/IntType$(Rl)-$(Ru)$(syn)$(Ni)$(en).png")
    plot(title=tl,xlabel="Strength of interactions",ylabel="Number of ecosystems")
    histogram!(log10.(mean1),fillalpha=0.75,label="Competiton")
    histogram!(log10.(mean2),fillalpha=0.75,label="Facilitation")
    histogram!(log10.(mean3),fillalpha=0.75,label="Syntrophy")
    histogram!(log10.(mean4),fillalpha=0.75,label="Pollution")
    savefig("Output/Fig5/IntStrength$(Rl)-$(Ru)$(syn)$(Ni)$(en).png")
    # Make range of ticks to label
    rgn = collect(-14:2:-2)
    ergn = fill("",length(rgn))
    # Find corresponding exponentials
    for i = 1:length(rgn)
        ergn[i] = L"10^{%$(rgn[i])}"
    end
    plot(title=tl,xlabel="Interaction strength",ylabel="Number of interactions")
    plot!(xticks = (rgn, ergn))
    histogram!(log10.(sts1),fillalpha=0.75,label="Competiton")
    histogram!(log10.(sts2),fillalpha=0.75,label="Facilitation")
    histogram!(log10.(sts3),fillalpha=0.75,label="Syntrophy")
    histogram!(log10.(sts4),fillalpha=0.75,label="Pollution")
    savefig("Output/Fig5/AllIntStrength$(Rl)-$(Ru)$(syn)$(Ni)$(en).png")
    return(nothing)
end

@time figure5(1,5,true,250,250,"i")
