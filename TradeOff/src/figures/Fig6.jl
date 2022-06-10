# Script to plot figure 6
using TradeOff
using Plots
using JLD
using ColorSchemes
using Plots.PlotMeasures
using LaTeXStrings
import PyPlot

function find_label(sim_type::Int64)
    # Assign label based on simulation type
    if sim_type == 1
        lb = "high free-energy"
    else
        lb = "low free-energy"
    end
    return(lb)
end

function figure6(rps::Int64,ims::Int64)
    println("Compiled")
    # Initialise plotting
    pyplot(dpi=200)
    # Load in color schemes
    a = ColorSchemes.Dark2_4.colors
    # Extract specific 4 colors from a color scheme
    bt = ColorSchemes.tab10.colors
    b = [bt[10];bt[9];bt[7];bt[5]]
    # Make plot objects
    p1 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),legend=:right,ylabel=L"\eta",title="ATP yield")
    p2 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),ylim=(0.55,Inf),legend=false)
    plot!(p2,ylabel="Fraction of free-energy transduced",title="Average reaction efficiency")
    p3 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),ylim=(0.25,Inf),legend=:bottomright)
    plot!(p3,ylabel="Fraction of free-energy transduced",title="Reaction efficency by reaction type")
    p4 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),ylim=(2.0,Inf),legend=false)
    plot!(p4,ylabel="Average number of reaction steps",title="Relative frequency of reaction types")
    # Loop over the 4 conditions
    for i = 1:2
        # Extract other simulation parameters from the function
        Np, Nt, M, d, μrange = sim_paras(i)
        # Read in appropriate files
        pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Paras$(ims)Ims.jld"
        if ~isfile(pfile)
            error("$(ims) immigrations run $(rN) is missing a parameter file")
        end
        # Find file name to load in
        tfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/RunStats$(ims)Ims.jld"
        # Check it actually exists
        if ~isfile(tfile)
            error("missing stats file for $(ims) immigrations simulations")
        end
        # Read in relevant data
        ps = load(pfile,"ps")
        # Now load out the times, and number of trajectories
        times = load(tfile,"times")
        no_via = load(tfile,"no_via")
        # Load in averages
        mn_via_η = load(tfile,"mn_via_η")
        mn_fr_ΔG = load(tfile,"mn_fr_ΔG")
        mn_fr_ΔG1 = load(tfile,"mn_fr_ΔG1")
        mn_fr_ΔG2 = load(tfile,"mn_fr_ΔG2")
        mn_av_steps = load(tfile,"mn_av_steps")
        # Load in standard deviations
        sd_via_η = load(tfile,"sd_via_η")
        sd_fr_ΔG = load(tfile,"sd_fr_ΔG")
        sd_fr_ΔG1 = load(tfile,"sd_fr_ΔG1")
        sd_fr_ΔG2 = load(tfile,"sd_fr_ΔG2")
        sd_av_steps = load(tfile,"sd_av_steps")
        # Calculate relevant standard errors
        se_via_η = sd_via_η./sqrt.(no_via)
        se_fr_ΔG = sd_fr_ΔG./sqrt.(no_via)
        se_fr_ΔG1 = sd_fr_ΔG1./sqrt.(no_via)
        se_fr_ΔG2 = sd_fr_ΔG2./sqrt.(no_via)
        se_av_steps = sd_av_steps./sqrt.(no_via)
        # Find appropriate label
        lb = find_label(i)
        # Plot the data to the relevant plot objects
        plot!(p1,times,mn_via_η,ribbon=se_via_η,color=a[i],label=lb)
        plot!(p2,times,mn_fr_ΔG,ribbon=se_fr_ΔG,color=a[i])
        plot!(p3,times,mn_fr_ΔG1,ribbon=se_fr_ΔG1,color=b[i],label="$(lb) 1-step")
        plot!(p3,times,mn_fr_ΔG2,ribbon=se_fr_ΔG2,color=b[2+i],label="$(lb) 2-step")
        plot!(p4,times,mn_av_steps,ribbon=se_av_steps,color=a[i])
    end
    # Check if directory exists and if not make it
    if ~isdir("Output/Fig6")
        mkdir("Output/Fig6")
    end
    # Save figures to this directory
    savefig(p1,"Output/Fig6/Eta.png")
    savefig(p2,"Output/Fig6/Efficency.png")
    savefig(p3,"Output/Fig6/CompEfficencies.png")
    savefig(p4,"Output/Fig6/Steps.png")
    # Add annotations
    px, py = annpos([0.0;5e7],[0.0;6.0],0.075,0.0)
    annotate!(p1,px,py,text("A",17,:black))
    px, py = annpos([0.0;5e7],[0.55;0.9],0.075,0.0)
    annotate!(p2,px,py,text("B",17,:black))
    px, py = annpos([0.0;5e7],[0.0;0.85],0.075,0.0)
    annotate!(p3,px,py,text("C",17,:black))
    px, py = annpos([0.0;5e7],[2.0;3.25],0.075,0.0)
    annotate!(p4,px,py,text("D",17,:black))
    # Plot all graphs as a single figure
    pt = plot(p1,p2,p3,p4,layout=(2,2),size=(1200,800),margin=5.0mm)
    savefig(pt,"Output/Fig6/figure6.png")
    return(nothing)
end

@time figure6(250,500)
