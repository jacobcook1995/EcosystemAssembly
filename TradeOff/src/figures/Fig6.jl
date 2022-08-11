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
    p5 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),ylim=(-Inf,Inf),legend=:bottomright)
    plot!(p5,ylabel="Fraction of free-energy transduced",title="Test plot")
    p6 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),ylim=(0.0,1.2),legend=:bottomright)
    plot!(p6,ylabel="Test plot",title="Test plot")
    p7 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),ylim=(0.0,1.2),legend=:bottomright)
    plot!(p7,ylabel="Test plot",title="Test plot")
    p8 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),ylim=(0.0,Inf),legend=:bottomright)
    plot!(p8,ylabel="Test plot",title="Test plot")
    # Loop over the 4 conditions
    for i = 3#1:2
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
        mn_η_stp = load(tfile,"mn_η_stp")
        mn_fr_ΔG = load(tfile,"mn_fr_ΔG")
        mn_fr_ΔG_stp = load(tfile,"mn_fr_ΔG_stp")
        mn_ϕP_stp = load(tfile,"mn_ϕP_stp")
        mn_av_steps = load(tfile,"mn_av_steps")
        # Load in standard deviations
        sd_via_η = load(tfile,"sd_via_η")
        sd_η_stp = load(tfile,"sd_η_stp")
        sd_fr_ΔG = load(tfile,"sd_fr_ΔG")
        sd_fr_ΔG_stp = load(tfile,"sd_fr_ΔG_stp")
        sd_ϕP_stp = load(tfile,"sd_ϕP_stp")
        sd_av_steps = load(tfile,"sd_av_steps")
        # Calculate relevant standard errors
        se_via_η = sd_via_η./sqrt.(no_via)
        se_η_stp = sd_η_stp./sqrt.(no_via)
        se_fr_ΔG = sd_fr_ΔG./sqrt.(no_via)
        se_fr_ΔG_stp  = sd_fr_ΔG_stp ./sqrt.(no_via)
        se_ϕP_stp  = sd_ϕP_stp ./sqrt.(no_via)
        se_av_steps = sd_av_steps./sqrt.(no_via)
        # Calculated expected eta value
        η_exp = mn_ϕP_stp[:,1].*mn_η_stp[:,1] .+ mn_ϕP_stp[:,2].*mn_η_stp[:,2]
        η_exp = η_exp .+ mn_ϕP_stp[:,3].*mn_η_stp[:,3] .+ mn_ϕP_stp[:,4].*mn_η_stp[:,4]
        # Find appropriate label
        lb = find_label(i)
        # Plot the data to the relevant plot objects
        plot!(p1,times,mn_via_η,ribbon=se_via_η,color=a[i],label=lb)
        plot!(p2,times,mn_fr_ΔG,ribbon=se_fr_ΔG,color=a[i])
        plot!(p3,times,mn_fr_ΔG_stp[:,1],ribbon=se_fr_ΔG_stp,color=b[1],label="$(lb) 1-step")
        plot!(p3,times,mn_fr_ΔG_stp[:,2],ribbon=se_fr_ΔG_stp,color=b[2],label="$(lb) 2-step")
        plot!(p3,times,mn_fr_ΔG_stp[:,3],ribbon=se_fr_ΔG_stp,color=b[3],label="$(lb) 3-step")
        plot!(p3,times,mn_fr_ΔG_stp[:,4],ribbon=se_fr_ΔG_stp,color=b[4],label="$(lb) 4-step")
        plot!(p4,times,mn_av_steps,ribbon=se_av_steps,color=a[i])
        plot!(p5,times,mn_η_stp[:,1],ribbon=se_η_stp,color=b[1],label="$(lb) 1-step")
        plot!(p5,times,mn_η_stp[:,2],ribbon=se_η_stp,color=b[2],label="$(lb) 2-step")
        plot!(p5,times,mn_η_stp[:,3],ribbon=se_η_stp,color=b[3],label="$(lb) 3-step")
        plot!(p5,times,mn_η_stp[:,4],ribbon=se_η_stp,color=b[4],label="$(lb) 4-step")
        plot!(p6,times,sum(mn_ϕP_stp[:,:],dims=2))
        plot!(p7,times,mn_ϕP_stp[:,1],ribbon=se_ϕP_stp,color=b[1],label="$(lb) 1-step")
        plot!(p7,times,mn_ϕP_stp[:,2],ribbon=se_ϕP_stp,color=b[2],label="$(lb) 2-step")
        plot!(p7,times,mn_ϕP_stp[:,3],ribbon=se_ϕP_stp,color=b[3],label="$(lb) 3-step")
        plot!(p7,times,mn_ϕP_stp[:,4],ribbon=se_ϕP_stp,color=b[4],label="$(lb) 4-step")
        plot!(p8,times,η_exp,label="")
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
    savefig(p5,"Output/Fig6/EtabyStep.png")
    savefig(p6,"Output/Fig6/Testplot.png")
    savefig(p7,"Output/Fig6/Comm_phi.png")
    savefig(p8,"Output/Fig6/Eta_comp.png")
    # Add annotations
    px, py = annpos([0.0;5e7],[0.0;6.0],0.075,0.0)
    annotate!(p1,[px,py,text("A",17,:black)])
    px, py = annpos([0.0;5e7],[0.55;0.9],0.075,0.0)
    annotate!(p2,[px,py,text("B",17,:black)])
    px, py = annpos([0.0;5e7],[0.0;0.85],0.075,0.0)
    annotate!(p3,[px,py,text("C",17,:black)])
    px, py = annpos([0.0;5e7],[2.0;3.25],0.075,0.0)
    annotate!(p4,[px,py,text("D",17,:black)])
    # Plot all graphs as a single figure
    pt = plot(p1,p2,p3,p4,layout=(2,2),size=(1200,800),margin=5.0mm)
    savefig(pt,"Output/Fig6/figure6.png")
    return(nothing)
end

@time figure6(250,500)
