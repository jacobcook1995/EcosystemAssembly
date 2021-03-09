# Script to construct figure 5
using Assembly
using Plots
using JLD
using LaTeXStrings
using Statistics
import PyPlot

function figure5(Rl::Int64,Ru::Int64,syns::Array{Bool,1},rps::Int64,Ni::Int64,en::String)
    println("Compiled!")
    # Preallocate memory to store number of interactions
    ins1 = zeros(rps,length(syns))
    ins2 = zeros(rps,length(syns))
    ins3 = zeros(rps,length(syns))
    ins4 = zeros(rps,length(syns))
    # Preallocate memory to store mean interaction strengths
    mean1 = fill(NaN,rps,length(syns))
    mean2 = fill(NaN,rps,length(syns))
    mean3 = fill(NaN,rps,length(syns))
    mean4 = fill(NaN,rps,length(syns))
    # Preallocate vectors to store all interaction strengths
    sts1 = Array{Array{Float64,1},1}(undef,2)
    sts2 = Array{Array{Float64,1},1}(undef,2)
    sts3 = Array{Array{Float64,1},1}(undef,2)
    sts4 = Array{Array{Float64,1},1}(undef,2)
    # Loop over parameter sets
    for i = 1:rps
        for j = 1:length(syns)
            # Read in relevant files
            ifile = "Data/$(Rl)-$(Ru)$(syns[j])$(Ni)$(en)/IntsReacs$(Rl)-$(Ru)Syn$(syns[j])Run$(i)Ns$(Ni).jld"
            if ~isfile(ifile)
                error("run $(i) is missing an interaction file")
            end
            # Just need to load ints and their strengths
            ints = load(ifile,"ints")
            in_str = load(ifile,"in_str")
            # Store number of each interaction type
            ins1[i,j] = count(x->x==1,ints)
            ins2[i,j] = count(x->x==2,ints)
            ins3[i,j] = count(x->x==3,ints)
            ins4[i,j] = count(x->x==4,ints)
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
                mean1[i,j] = mean(str1)
                # If vector doesn't already exist create it
                if ~isdefined(sts1,j)
                    sts1[j] = str1
                else
                    sts1[j] = cat(sts1[j],str1,dims=1)
                end
            end
            # Do for each interaction type
            if length(str2) >= 1
                mean2[i,j] = mean(str2)
                # If vector doesn't already exist create it
                if ~isdefined(sts2,j)
                    sts2[j] = str2
                else
                    sts2[j] = cat(sts2[j],str2,dims=1)
                end
            end
            if length(str3) >= 1
                mean3[i,j] = mean(str3)
                # If vector doesn't already exist create it
                if ~isdefined(sts3,j)
                    sts3[j] = str3
                else
                    sts3[j] = cat(sts3[j],str3,dims=1)
                end
            end
            if length(str4) >= 1
                mean4[i,j] = mean(str4)
                # If vector doesn't already exist create it
                if ~isdefined(sts4,j)
                    sts4[j] = str4
                else
                    sts4[j] = cat(sts4[j],str4,dims=1)
                end
            end
        end
    end
    # Setup plotting
    pyplot()
    theme(:wong2,dpi=200)
    # Loop over syntrophy conditions
    for j = 1:length(syns)
        # Make plot title
        tl = ""
        if syns[j] == true
            tl = "$(Rl)-$(Ru) reactions per strain"
        else
            tl = "$(Rl)-$(Ru) reactions per strain (no syntrophy)"
        end
        # Now plot both interactions types and strengths
        plot(title=tl,xlabel="Number of interactions",ylabel="Number of ecosystems")
        histogram!(ins1[:,j],fillalpha=0.75,label="Competiton")
        histogram!(ins2[:,j],fillalpha=0.75,label="Facilitation")
        histogram!(ins3[:,j],fillalpha=0.75,label="Syntrophy")
        histogram!(ins4[:,j],fillalpha=0.75,label="Pollution")
        savefig("Output/Fig5/IntType$(Rl)-$(Ru)$(syns[j])$(Ni)$(en).png")
        plot(title=tl,xlabel="Strength of interactions",ylabel="Number of ecosystems")
        histogram!(log10.(mean1[:,j]),fillalpha=0.75,label="Competiton")
        histogram!(log10.(mean2[:,j]),fillalpha=0.75,label="Facilitation")
        histogram!(log10.(mean3[:,j]),fillalpha=0.75,label="Syntrophy")
        histogram!(log10.(mean4[:,j]),fillalpha=0.75,label="Pollution")
        savefig("Output/Fig5/IntStrength$(Rl)-$(Ru)$(syns[j])$(Ni)$(en).png")
        # Make range of ticks to label
        rgn = collect(-14:2:-2)
        ergn = fill("",length(rgn))
        # Find corresponding exponentials
        for i = 1:length(rgn)
            ergn[i] = L"10^{%$(rgn[i])}"
        end
        plot(title=tl,xlabel="Interaction strength",ylabel="Number of interactions")
        plot!(xticks = (rgn, ergn))
        histogram!(log10.(sts1[j]),fillalpha=0.75,label="Competiton")
        histogram!(log10.(sts2[j]),fillalpha=0.75,label="Facilitation")
        histogram!(log10.(sts3[j]),fillalpha=0.75,label="Syntrophy")
        histogram!(log10.(sts4[j]),fillalpha=0.75,label="Pollution")
        savefig("Output/Fig5/AllIntStrength$(Rl)-$(Ru)$(syns[j])$(Ni)$(en).png")
    end
    return(nothing)
end

@time figure5(1,5,[true,false],250,250,"i")
