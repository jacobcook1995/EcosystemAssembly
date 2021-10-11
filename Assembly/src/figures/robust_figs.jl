# Script to make robustness figures
using Assembly
using Plots
using JLD
using StatsPlots
using StatsBase
using Statistics
using Plots.PlotMeasures
using LaTeXStrings
using ColorSchemes
using KernelDensity
import PyPlot

# function to find average efficency
function av_eff(pop::Array{Float64,1},conc::Array{Float64,1},ms::Array{MicrobeP,1},ps::FullParameters)
    # Define mimimum product to substrate ratio (to calculate) the efficency
    mr = 1e-2
    # Find indices of surving populations
    inds = findall(x->x>0.0,pop)
    # Set initial efficency value
    efT = 0.0
    # Loop over indices
    for i = inds
        # Loop over number of reactions
        for j = 1:ms[i].R
            # Calculate proportional effiency
            eff = -(ms[i].η[j]*ΔGATP)/(ps.reacs[ms[i].Reacs[j]].ΔG0+Rgas*ps.T*log(mr))
            # Add this to total weighted by population, and metabolic protein fraction
            efT += pop[i]*ms[i].ϕP[j]*eff
        end
    end
    # Average across populations
    avef = efT/sum(pop)
    return(avef)
end

# funtion to find average growth rate
function av_λ(pop::Array{Float64,1},as::Array{Float64,1},ϕRs::Array{Float64,1},ms::Array{MicrobeP,1})
    # Find indices of surving populations
    inds = findall(x->x>0.0,pop)
    # weighted total growth rate (starts at zero)
    λT = 0.0
    # Loop over survivors
    for i = inds
        # Find growth rate for this particular species
        λt = λs(as[i],ϕRs[i],ms[i])
        # Weight by population and add to total
        λT += pop[i]*λt
    end
    # Divide weighted growth rate by total abundance
    λT /= sum(pop)
    return(λT)
end

function rob_figs(Rls::Array{Int64,1},Rus::Array{Int64,1},syns::Array{Bool,1},ens::Array{String,1},
                Ni::Int64,Nr::Int64,Nr_1::Int64)
    # Check if all these vectors are the same length
    if length(Rls) != length(Rus) || length(Rls) != length(syns) || length(Rls) != length(ens)
        error("length of vectors doesn't match")
    end
    println("Compiled!")
    # Count number of parameter sets (for each condition)
    Ns = length(Rls)
    # Containers to store number of survivors
    svs = zeros(Int64,Ns,Nr)
    svs_1 = zeros(Int64,Ns,Nr_1)
    # Setup plotting
    pyplot()
    theme(:wong2,dpi=300)
    wongc = wong2_palette()
    # Find corresponding indices of reordered labels
    pos = zeros(Float64,Ns)
    c = Array{RGBA,1}(undef,Ns)
    # Find posistions for mean and std
    for i = 1:Ns
        p = 1.0
        if syns[i] == true
            p += 0.5
            # Set color
            c[i] = wongc[1]
        else
            # Choose nice grey as "light" color
            c[i] = RGB(([170, 170, 170] / 255)...)
        end
        if ens[i] == "l"
            p += 1.5
        end
        pos[i] = p
    end
    # Set titles for the six new conditions
    ttls = fill("",6)
    ttls[1] = "Increase $(L"n_P")"
    ttls[2] = "Increase $(L"f_b")"
    ttls[3] = "Increase $(L"\phi_Q")"
    ttls[4] = "Increase $(L"\gamma_{\frac{1}{2}}")"
    ttls[5] = "High saturation"
    ttls[6] = "Low saturation"
    # Set whether the six new conditions show significant differences (aside the obvious one)
    sig_dif = fill(false,6)
    sig_dif[1] = true
    sig_dif[3] = true
    sig_dif[5] = true
    # Load in data for the first figure by looping over parameter sets
    for i = 1:Ns
        for j = 1:Nr_1
            # Read in relevant files
            pfile = "Data/$(Rls[i])-$(Rus[i])$(syns[i])$(Ni)$(ens[i])/RedParasReacs$(Rls[i])-$(Rus[i])Syn$(syns[i])Run$(j)Ns$(Ni).jld"
            if ~isfile(pfile)
                error("parameter set $(i) run $(j) is missing a parameter file")
            end
            ofile = "Data/$(Rls[i])-$(Rus[i])$(syns[i])$(Ni)$(ens[i])/RedOutputReacs$(Rls[i])-$(Rus[i])Syn$(syns[i])Run$(j)Ns$(Ni).jld"
            if ~isfile(ofile)
                error("parameter set $(i) run $(j) is missing an output file")
            end
            # Only want final parameter set
            ps = load(pfile,"ps")
            inf_out = load(ofile,"inf_out")
            # Save number of survivors
            svs_1[i,j] = ps.N
        end
    end
    # Make first figure
    p1 = plot(ylabel="Number of surviving species",xlim=(0.5,3.5),xlabel="Energy supply")
    plot!(p1,xticks=([1.25,2.75],["high","low"]),title="Original case",legend=:right)
    # Plot means
    for i = 1:Ns
        # Calculate mean
        mn = mean(svs_1[i,:])
        # Calculate 99% confidence interval
        sdn = sem(svs_1[i,:])*2.576
        # Label empty
        lb = ""
        # Unless energy is high
        if ens[i] == "h"
            if syns[i] == true
                lb = "Reversible"
            else
                lb = "M–M"
            end
        end
        scatter!(p1,[pos[i]],[mn],yerror=[sdn],label=lb,color=c[i],ms=6,msc=c[i])
    end
    # Add bracket for significance plot
    plot!(p1,[2.5,3.0],[5.0,5.0],linecolor=:black,label="")
    plot!(p1,[2.5,2.5],[4.6,5.01],linecolor=:black,label="")
    plot!(p1,[3.0,3.0],[4.6,5.01],linecolor=:black,label="")
    # Then add star above the bracket
    scatter!(p1,[2.75],[5.25],color=:black,shape=:star6,label="")
    # Split into two plots by changing limits
    p1a = plot!(p1,ylim=(0.0,14.0))
    p1b = plot!(deepcopy(p1),ylim=(2.0,17.0))
    # Preallocate vector of subplots
    p = Array{Plots.Plot,1}(undef,6)
    # Assign basic plot features for each subplot
    for i = 1:6
        p[i] = plot(xlim=(0.5,3.5),xlabel="Energy supply")
    end
    # Loop over the 6 different conditions
    for l = 1:6
        # Find title of option
        title_op = options_titles(l)
        # Loop over parameter sets
        for i = 1:Ns
            for j = 1:Nr
                # Read in relevant files
                pfile = "Paras/$(title_op)/$(Rls[i])-$(Rus[i])$(syns[i])$(ens[i])/RedParasReacs$(Rls[i])-$(Rus[i])Syn$(syns[i])Run$(j)Ns$(Ni).jld"
                if ~isfile(pfile)
                    error("parameter set $(i) run $(j) is missing a parameter file")
                end
                ofile = "Data/$(title_op)/$(Rls[i])-$(Rus[i])$(syns[i])$(ens[i])/RedOutputReacs$(Rls[i])-$(Rus[i])Syn$(syns[i])Run$(j)Ns$(Ni).jld"
                if ~isfile(ofile)
                    error("parameter set $(i) run $(j) is missing an output file")
                end
                # Only want final parameter set
                ps = load(pfile,"ps")
                inf_out = load(ofile,"inf_out")
                # Save number of survivors
                svs[i,j] = ps.N
            end
        end
        # Container to store mean + sd for each case
        msd = zeros(Ns)
        plot!(p[l],xticks=([1.25,2.75],["high","low"]),title=ttls[l])
        # Plot means
        for i = 1:Ns
            # Calculate mean
            mn = mean(svs[i,:])
            # Calculate 99% confidence interval
            sdn = sem(svs[i,:])*2.576
            scatter!(p[l],[pos[i]],[mn],yerror=[sdn],label="",color=c[i],ms=6,msc=c[i])
        end
        # Check if there's a significant difference between the two low free-energy conditions
        if sig_dif[l] == true
            # If so add bracket for significance plot
            plot!(p[l],[2.5,3.0],[5.0,5.0],linecolor=:black,label="")
            plot!(p[l],[2.5,2.5],[4.6,5.01],linecolor=:black,label="")
            plot!(p[l],[3.0,3.0],[4.6,5.01],linecolor=:black,label="")
            # Then add star above the bracket
            scatter!(p[l],[2.75],[5.25],color=:black,shape=:star6,label="")
        end
        # Now set ylims
        if l <= 3
            p[l] = plot!(p[l],ylim=(0.0,14.0))
        else
            p[l] = plot!(p[l],ylim=(2.0,17.0))
        end
    end
    # Combine plots so that comparions can be performed
    pc1 = plot(p1a,p[1],p[2],p[3],layout=(1,4),size=(800,400),guidefontsize=13,legendfontsize=8,tickfontsize=11)
    savefig(pc1,"Output/SI/Surv_comp_1.png")
    pc2 = plot(p1b,p[4],p[5],p[6],layout=(1,4),size=(800,400),guidefontsize=13,legendfontsize=8,tickfontsize=11)
    savefig(pc2,"Output/SI/Surv_comp_2.png")
    return(nothing)
end

# Hard code parameters here
l = [1,1,1,1]
u = [5,5,5,5]
s = [true,true,false,false]
e = ["l","h","l","h"]

@time rob_figs(l,u,s,e,250,50,250)
