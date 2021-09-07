# Script to plot figure 6
using TradeOff
using Plots
using JLD
using ColorSchemes
using Plots.PlotMeasures
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

function figure6(rps::Int64,ims::Int64)
    println("Compiled")
    # Initialise plotting
    pyplot(dpi=200)
    # Load in color scheme
    a = ColorSchemes.Dark2_4.colors
    # Make plot objects
    p1 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),legend=:bottomright,ylabel="Number of strains")
    p2 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),legend=false,ylabel="Number of viable strains")
    p3 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),legend=false,ylabel="Proportion of strains stable")
    p4 = plot(xlabel="Times (s)",xlim=(-Inf,5e7),legend=false,ylabel="Proportion of invaders that are feasible",ylims=(-Inf,0.5))
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
        sfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/SnapDataStats$(ims)Ims.jld"
        if ~isfile(sfile)
            error("$(ims) immigrations run $(rN) is missing an output file")
        end
        # Read in relevant data
        ps = load(pfile,"ps")
        # Now load out the times, and number of trajectories
        times_t = load(tfile,"times")
        no_sims = load(tfile,"no_sims")
        no_via = load(tfile,"no_via")
        # Load in averages
        mn_svt = load(tfile,"mn_svt")
        mn_tsvt = load(tfile,"mn_tsvt")
        # Load in standard deviations
        sd_svt = load(tfile,"sd_svt")
        sd_tsvt = load(tfile,"sd_tsvt")
        # Load in times fro snapshot data
        times_s = load(sfile,"times")
        # Load in all the rest of the relevant snapshot data
        gp = load(sfile,"gp")
        mn_stb = load(sfile,"mn_stb")
        mn_inc = load(sfile,"mn_inc")
        mn_dec = load(sfile,"mn_dec")
        sd_stb = load(sfile,"sd_stb")
        sd_inc = load(sfile,"sd_inc")
        sd_dec = load(sfile,"sd_dec")
        st_r = load(sfile,"st_r")
        # Calculate relevant standard errors
        se_svt = sd_svt./sqrt.(no_sims)
        se_tsvt = sd_tsvt./sqrt.(no_sims)
        se_stb = sd_stb./sqrt.(st_r)
        se_inc = sd_inc./sqrt.(st_r)
        se_dec = sd_dec./sqrt.(st_r)
        # Calculate proportion stable
        mn_props = mn_stb./(mn_stb.+mn_inc.+mn_dec)
        # Do error propogation
        se_props = sqrt.((se_stb.^2).*(mn_inc.+mn_dec).^2 + (se_dec.^2 .+ se_inc.^2).*mn_stb.^2)/((mn_stb+mn_dec+mn_inc).^2)
        # Find appropriate label
        lb = find_label(i)
        # Plot the data to the relevant plot objects
        plot!(p1,times_t,mn_svt,ribbon=se_svt,label=lb,color=a[i])
        plot!(p2,times_t,mn_tsvt,ribbon=se_tsvt,color=a[i])
        plot!(p3,times_s,mn_props,ribbon=se_props,color=a[i])
        plot!(p4,times_s,gp,color=a[i])
    end
    # Check if directory exists and if not make it
    if ~isdir("Output/Fig6")
        mkdir("Output/Fig6")
    end
    # Save figures to this directory
    savefig(p1,"Output/Fig6/TotalStrains.png")
    savefig(p2,"Output/Fig6/ViableStrains.png")
    savefig(p3,"Output/Fig6/Stable.png")
    savefig(p4,"Output/Fig6/Feasible.png")
    # Add annotations
    px, py = annpos([0.0;5e7],[0.0;35.0],0.075,0.0)
    annotate!(p1,px,py,text("A",17,:black))
    px, py = annpos([0.0;5e7],[0.0;22.0],0.075,0.0)
    annotate!(p2,px,py,text("B",17,:black))
    px, py = annpos([0.0;5e7],[0.0;0.82],0.075,0.0)
    annotate!(p3,px,py,text("C",17,:black))
    px, py = annpos([0.0;5e7],[0.0;0.5],0.075,0.0)
    annotate!(p4,px,py,text("D",17,:black))
    # Plot all graphs as a single figure
    pt = plot(p1,p2,p3,p4,layout=4,size=(1200,800),margin=5.0mm)
    savefig(pt,"Output/Fig6/figure6.png")
    return(nothing)
end

@time figure6(250,500)
