# Script to plot figure 7
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
        lb = "high free energy + low loss"
    elseif sim_type == 2
        lb = "low free energy + low loss"
    elseif sim_type == 3
        lb = "high free energy + high loss"
    else
        lb = "low free energy + high loss"
    end
    return(lb)
end

function figure8(rps::Int64,ims::Int64)
    println("Compiled")
    # Initialise plotting
    pyplot(dpi=200)
    # Load in color scheme
    a = ColorSchemes.Dark2_4.colors
    # Make plot objects
    p1 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),legend=:topright,ylabel=L"\omega")
    p2 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),legend=false,ylabel=L"\phi_R")
    p3 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),legend=false,ylabel=L"\eta")
    p4 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),legend=false,ylabel="Fraction transduced")
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
        mn_via_ω = load(tfile,"mn_via_ω")
        mn_via_ϕR = load(tfile,"mn_via_ϕR")
        mn_via_η = load(tfile,"mn_via_η")
        mn_fr_ΔG = load(tfile,"mn_fr_ΔG")
        # Load in standard deviations
        sd_via_ω = load(tfile,"sd_via_ω")
        sd_via_ϕR = load(tfile,"sd_via_ϕR")
        sd_via_η = load(tfile,"sd_via_η")
        sd_fr_ΔG = load(tfile,"sd_fr_ΔG")
        # Calculate relevant standard errors
        se_via_ω = sd_via_ω./sqrt.(no_via)
        se_via_ϕR = sd_via_ϕR./sqrt.(no_via)
        se_via_η = sd_via_η./sqrt.(no_via)
        se_fr_ΔG = sd_fr_ΔG./sqrt.(no_via)
        # Find appropriate label
        lb = find_label(i)
        # Plot the data to the relevant plot objects
        plot!(p1,times,mn_via_ω,ribbon=se_via_ω,label=lb,color=a[i])
        plot!(p2,times,mn_via_ϕR,ribbon=se_via_ϕR,color=a[i])
        plot!(p3,times,mn_via_η,ribbon=se_via_η,color=a[i])
        plot!(p4,times,mn_fr_ΔG,ribbon=se_fr_ΔG,color=a[i])
    end
    # Check if directory exists and if not make it
    if ~isdir("Output/Fig8")
        mkdir("Output/Fig8")
    end
    # Save figures to this directory
    savefig(p1,"Output/Fig8/omegas.png")
    savefig(p2,"Output/Fig8/Rfracs.png")
    savefig(p3,"Output/Fig8/Eta.png")
    savefig(p4,"Output/Fig8/Efficency.png")
    # Add annotations
    px, py = annpos([0.0;5e7],[0.0,1.0],0.075,0.0)
    annotate!(p1,px,py,text("A",17,:black))
    # 32 over shoots
    px, py = annpos([0.0;5e7],[0.0;0.305],0.075,0.0)
    annotate!(p2,px,py,text("B",17,:black))
    px, py = annpos([0.0;5e7],[0.0;5.75],0.075,0.0)
    annotate!(p3,px,py,text("C",17,:black))
    # 1.0 is too low
    px, py = annpos([0.0;5e7],[0.0;1.04],0.075,0.0)
    annotate!(p4,px,py,text("D",17,:black))
    # Plot all graphs as a single figure
    pt = plot(p1,p2,p3,p4,layout=(2,2),size=(1200,800),margin=5.0mm)
    savefig(pt,"Output/Fig8/figure8.png")
    return(nothing)
end

@time figure8(250,500)
