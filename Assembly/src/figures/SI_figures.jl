# Script to construct figure 5
using Assembly
using Plots
using JLD
using LaTeXStrings
using StatsBase
using Plots.PlotMeasures
import PyPlot

# Make figure plots for the interactions
function SI_ints(Rl::Int64,Ru::Int64,syns::Array{Bool,1},rps::Int64,Ni::Int64,en::String)
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
    # Make array of bools to check if data has already been saved
    fll = fill(false,4,2)
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
                if fll[1,j] == false
                    sts1[j] = str1
                    # Update counter
                    fll[1,j] = true
                else
                    sts1[j] = cat(sts1[j],str1,dims=1)
                end
            end
            # Do for each interaction type
            if length(str2) >= 1
                mean2[i,j] = mean(str2)
                # If vector doesn't already exist create it
                if fll[2,j] == false
                    sts2[j] = str2
                    # Update counter
                    fll[2,j] = true
                else
                    sts2[j] = cat(sts2[j],str2,dims=1)
                end
            end
            if length(str3) >= 1
                mean3[i,j] = mean(str3)
                # If vector doesn't already exist create it
                if fll[3,j] == false
                    sts3[j] = str3
                    # Update counter
                    fll[3,j] = true
                else
                    sts3[j] = cat(sts3[j],str3,dims=1)
                end
            end
            if length(str4) >= 1
                mean4[i,j] = mean(str4)
                # If vector doesn't already exist create it
                if fll[4,j] == false
                    sts4[j] = str4
                    # Update counter
                    fll[4,j] = true
                else
                    sts4[j] = cat(sts4[j],str4,dims=1)
                end
            end
        end
    end
    # Find percentage of each interactions type
    cvsf = 100.0*(ins1)./(ins1.+ins2)
    pvsn = 100.0*(ins2.+ins3)./(ins1.+ins2.+ins3.+ins4)
    tvnt = 100.0*(ins3.+ins4)./(ins1.+ins2.+ins3.+ins4)
    # Setup plotting
    pyplot()
    theme(:wong2,dpi=300)
    # Make a color scheme
    pl = Array{RGBA,1}(undef,3)
    # Define the 4 colors I want to show
    pl[1] = RGB(([254,27,182] / 255)...) # red
    pl[2] = RGB(([182,254,27] / 255)...) # green
    pl[3] = RGB(([25,181,254] / 255)...) # blue
    # Preallocate array to store plots
    p = Array{Plots.Plot,2}(undef,2,length(syns))
    # Loop over syntrophy conditions
    for j = 1:length(syns)
        # Make plot title
        tl = ""
        if syns[j] == true
            tl = "Reversible kinetics"
        else
            tl = "Michaelis–Menten kinetics"
        end
        # Make bins for proportions
        pbins = range(-1.0,stop=101.0,length=52)
        # Now plot both interactions types and strengths
        p[1,j] = plot(title=tl,ylabel="Number of ecosystems")
        # Turn off legend for reversible case
        if syns[j] == true
            plot!(p[1,j],legend=false,xlabel="Percentage of interactions")
        end
        histogram!(p[1,j],cvsf[:,j],fillalpha=0.75,label="Competative",bins=pbins,color=pl[1])
        histogram!(p[1,j],pvsn[:,j],fillalpha=0.75,label="Positive",bins=pbins,color=pl[2])
        histogram!(p[1,j],tvnt[:,j],fillalpha=0.75,label="Thermodynamic",bins=pbins,color=pl[3])
        # Choose which letter to annotate
        if syns[j] == false
            # Then check energy supply
            if en == "h"
                # Add annotation
                px, py = annpos([0.0,100.0],[0.0,135.0],0.15,0.05)
                annotate!(p[1,j],px,py,text("A",17,:black))
            else
                # Add annotation
                px, py = annpos([0.0,100.0],[0.0,63.0],0.15,0.05)
                annotate!(p[1,j],px,py,text("A",17,:black))
            end
        else
            # Then check energy supply
            if en == "h"
                # Add annotation
                px, py = annpos([0.0,100.0],[0.0,67.5],0.15,0.05)
                annotate!(p[1,j],px,py,text("C",17,:black))
            else
                # Add annotation
                px, py = annpos([0.0,100.0],[0.0,42.0],0.15,0.05)
                annotate!(p[1,j],px,py,text("C",17,:black))
            end
        end
        savefig(p[1,j],"Output/SI/IntType$(Rl)-$(Ru)$(syns[j])$(Ni)$(en).png")
        # Make range of ticks to label
        rgn = collect(-15:3:0)
        ergn = fill("",length(rgn))
        # Find corresponding exponentials
        for i = 1:length(rgn)
            ergn[i] = L"10^{%$(rgn[i])}"
        end
        # Make bins for proportions
        sbins = range(-15.0,stop=0.0,length=52)
        # Now plot all strengths
        p[2,j] = plot(title=tl,ylabel="Number of interactions",legend=:topleft)
        plot!(p[2,j],xticks=(rgn,ergn))
        # Add xlabel for reversible case
        if syns[j] == true
            plot!(p[2,j],legend=false,xlabel="Interaction strength")
        end
        histogram!(p[2,j],log10.(sts1[j]),fillalpha=0.75,label="Competiton",bins=sbins)
        histogram!(p[2,j],log10.(sts2[j]),fillalpha=0.75,label="Facilitation",bins=sbins)
        histogram!(p[2,j],log10.(sts3[j]),fillalpha=0.75,label="Syntrophy",bins=sbins)
        histogram!(p[2,j],log10.(sts4[j]),fillalpha=0.75,label="Pollution",bins=sbins)
        # Choose which letter to annotate
        if syns[j] == false
            # Then check energy supply
            if en == "h"
                # Add annotation
                px, py = annpos([-15.0,0.0],[0.0,4300.0],0.15,0.05)
                annotate!(p[2,j],px,py,text("B",17,:black))
            else
                # Add annotation
                px, py = annpos([-15.0,0.0],[0.0,280.0],0.15,0.05)
                annotate!(p[2,j],px,py,text("B",17,:black))
            end
        else
            # Then check energy supply
            if en == "h"
                # Add annotation
                px, py = annpos([-15.0,0.0],[0.0,4100.0],0.15,0.05)
                annotate!(p[2,j],px,py,text("D",17,:black))
            else
                # Add annotation
                px, py = annpos([-15.0,0.0],[0.0,350.0],0.15,0.05)
                annotate!(p[2,j],px,py,text("D",17,:black))
            end

        end
        # Fit pollution histogram
        hp = fit(Histogram,log10.(sts4[j]),sbins,closed=:right)
        # Find maximum height of this pollution histogram
        hmaxp = maximum(hp.weights)
        # Find index of this maximum
        mind = findfirst(x->x==hmaxp,hp.weights)
        # Find first index after this that's less than half the value
        dind = findfirst(x->x<=0.5*hmaxp,hp.weights[(mind+1):end])
        # Extract bins from histogram
        rn = collect(sbins)
        # Next need to find x posistion of this half maximum point
        xpos = (rn[mind+dind]+rn[mind+dind-1])/2
        # Put in a vertical line here
        vline!(p[2,j],[xpos],color=:red,label="")
        savefig(p[2,j],"Output/SI/AllIntStrength$(Rl)-$(Ru)$(syns[j])$(Ni)$(en).png")
    end
    # Combine all graphs and save
    pt = plot(p[1,2],p[2,2],p[1,1],p[2,1],layout=4,size=(1200,800),margin=5.0mm)
    # Save as a different name based on energy level
    if en == "h"
        savefig(pt,"Output/SI/HighInts.png")
    elseif en == "l"
        savefig(pt,"Output/SI/LowInts.png")
    end
    return(nothing)
end
# NEED TO FIX LABELS!!!!!

# Run both high and low as want both plots for SI
@time SI_ints(1,5,[true,false],250,250,"h")
@time SI_ints(1,5,[true,false],250,250,"l")
