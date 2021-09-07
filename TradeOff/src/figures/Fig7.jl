# Script to plot figure 7
using TradeOff
using Plots
using JLD
using ColorSchemes
using Plots.PlotMeasures
import PyPlot

function figure7(rps::Int64,ims::Int64)
    println("Compiled")
    # Initialise plotting
    pyplot(dpi=200)
    # Load in color scheme
    a = ColorSchemes.tab20c.colors
    # Make plot objects
    p1 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),legend=:bottomright,ylabel="",title="Low biomass loss")
    p2 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),legend=false,ylabel="",title="High biomass loss")
    # Loop over the 4 conditions
    for i = 1:4
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
        mn_via_R = load(tfile,"mn_via_R")
        # Load in standard deviations
        sd_via_R = load(tfile,"sd_via_R")
        # Preallocate standard errors
        se_via_R = zeros(size(sd_via_R))
        # Calculation (slightly) different in the viable case
        for i = 1:size(sd_via_R,1)
            se_via_R[i,:] = sd_via_R[i,:]./sqrt.(no_via)
        end
        # Check if this is a low or high biomass loss case
        if i < 3
            # Check whether free energy is high or low
            if i == 1
                plot!(p1,times,mn_via_R[1,:],ribbon=se_via_R[1,:],label="R=1 high free-energy",color=a[1])
                plot!(p1,times,mn_via_R[7,:],ribbon=se_via_R[7,:],label="R=7 high free-energy",color=a[2])
            else
                plot!(p1,times,mn_via_R[1,:],ribbon=se_via_R[1,:],label="R=1 low free-energy",color=a[5])
                plot!(p1,times,mn_via_R[7,:],ribbon=se_via_R[7,:],label="R=7 low free-energy",color=a[6])
            end
        else
            # Check whether free energy is high or low
            if i == 3
                plot!(p2,times,mn_via_R[1,:],ribbon=se_via_R[1,:],color=a[1])
                plot!(p2,times,mn_via_R[7,:],ribbon=se_via_R[7,:],color=a[2])
            else
                plot!(p2,times,mn_via_R[1,:],ribbon=se_via_R[1,:],color=a[5])
                plot!(p2,times,mn_via_R[7,:],ribbon=se_via_R[7,:],color=a[6])
            end
        end
    end
    # Check if directory exists and if not make it
    if ~isdir("Output/Fig7")
        mkdir("Output/Fig7")
    end
    # Save figures to this directory
    savefig(p1,"Output/Fig7/LowLoss.png")
    savefig(p2,"Output/Fig7/HighLoss.png")
    # Add annotations
    px, py = annpos([0.0;5e7],[0.0;8.5],0.075,0.0)
    annotate!(p1,px,py,text("A",17,:black))
    px, py = annpos([0.0;5e7],[0.0;8.0],0.075,0.0)
    annotate!(p2,px,py,text("B",17,:black))
    # Now want to make a plot incorperating both previous plots
    pt = plot(p1,p2,layout=(2,1),size=(750,800),margin=5.0mm,rightmargin=20mm)
    savefig(pt,"Output/Fig7/figure7.png")
    return(nothing)
end

@time figure7(250,500)
