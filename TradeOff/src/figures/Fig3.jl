# Script to plot figure 3
using TradeOff
using Plots
using JLD
using Plots.PlotMeasures
import PyPlot

# STORE THIS HERE FOR NOW INCASE IT BECOMES USEFUL
function figure4(rps::Int64,ims::Int64,sim_type::Int64)
    println("Compiled")
    # Extract other simulation parameters from the function
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Read in appropriate files
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    sfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/SnapDataStats$(ims)Ims.jld"
    if ~isfile(sfile)
        error("$(ims) immigrations run $(rN) is missing an output file")
    end
    # Read in relevant data
    ps = load(pfile,"ps")
    # Save growth probabilities
    gp = load(sfile,"gp")
    mn_stb = load(sfile,"mn_stb")
    mn_inc = load(sfile,"mn_inc")
    mn_dec = load(sfile,"mn_dec")
    sd_stb = load(sfile,"sd_stb")
    sd_inc = load(sfile,"sd_inc")
    sd_dec = load(sfile,"sd_dec")
    st_r = load(sfile,"st_r")
    s_times = load(sfile,"times")
    # Find file name to load in
    tfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/RunStats$(ims)Ims.jld"
    # Check it actually exists
    if ~isfile(tfile)
        error("missing stats file for $(ims) immigrations simulations")
    end
    # Now load out the times, and number of trajectories
    t_times = load(tfile,"times")
    no_sims = load(tfile,"no_sims")
    no_via = load(tfile,"no_via")
    # Load in averages
    mn_svt = load(tfile,"mn_svt")
    mn_tsvt = load(tfile,"mn_tsvt")
    # Load in standard deviations
    sd_svt = load(tfile,"sd_svt")
    sd_tsvt = load(tfile,"sd_tsvt")
    println("Data read in")
    # Check if directory exists and if not make it
    if ~isdir("Output/Fig4")
        mkdir("Output/Fig4")
    end
    # Calculate standard errors
    se_stb = sd_stb./sqrt.(st_r)
    se_inc = sd_inc./sqrt.(st_r)
    se_dec = sd_dec./sqrt.(st_r)
    se_svt = sd_svt./sqrt.(no_sims)
    se_tsvt = sd_tsvt./sqrt.(no_sims)
    # Setup plotting
    pyplot(dpi=200)
    # Load in color scheme
    a = ColorSchemes.tab10.colors
    # Plot total survivors
    p1 = plot(xlabel="Time (s)",ylabel="Number of strains",xlim=(-Inf,5e7),legend=:topleft)
    plot!(p1,t_times,mn_svt,ribbon=se_svt,label="Total",color=a[1])
    plot!(p1,t_times,mn_tsvt,ribbon=se_tsvt,label="Viable",color=a[2])
    # Add annotation
    px, py = annpos([0.0;5e7],[0.0;35.0],0.05,0.0)
    annotate!(p1,px,py,text("A",17,:black))
    # Twin the x-axis
    pt = twinx(p1)
    # Ensure same limits are used
    plot!(xlim=(-Inf,5e7))
    # Then plot the proportion feasible
    p2 = plot!(pt,s_times,gp,label="Feasibility",color=a[3],legend=:topright,ylabel="Proportion of invaders that are feasible")
    savefig(p2,"Output/Fig4/Strains.png")
    # Then plot various means
    p3 = plot(xlabel="Time (s)",ylabel="Number of strains",xlim=(-Inf,5e7))
    plot!(p3,s_times,mn_stb,ribbon=se_stb,label="Stable",color=a[4])
    plot!(p3,s_times,mn_inc,ribbon=se_inc,label="Growing",color=a[5])
    plot!(p3,s_times,mn_dec,ribbon=se_dec,label="Declining",color=a[6])
    # Add annotation
    px, py = annpos([0.0;5e7],[0.0;21.0],0.05,0.0)
    annotate!(p3,px,py,text("B",17,:black))
    savefig(p3,"Output/Fig4/GrowvsDecl.png")
    # Now want to make a plot incorperating both previous plots
    pt = plot(p2,p3,layout=(2,1),size=(750,800),margin=5.0mm,rightmargin=20mm)
    savefig(pt,"Output/Fig4/figure4.png")
    return(nothing)
end

function figure3(rps::Int64,ims::Int64,sim_type::Int64)
    println("Compiled")
    # Extract other simulation parameters from the function
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Read in appropriate files
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    ofile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Run$(rN)Data$(ims)Ims.jld"
    if ~isfile(ofile)
        error("$(ims) immigrations run $(rN) is missing an output file")
    end
    # Read in relevant data
    ps = load(pfile,"ps")
    traj = load(ofile,"traj")
    T = load(ofile,"T")
    micd = load(ofile,"micd")
    its = load(ofile,"its")
    println("Data read in")
    # Find C from a function
    C = merge_data(ps,traj,T,micd,its)
    println("Data merged")
    # Check if directory exists and if not make it
    if ~isdir("Output/Fig3")
        mkdir("Output/Fig3")
    end
    # Find total number of strains
    totN = length(micd)
    # Find indices of extinct strains
    ext = isnan.(C[end,1:totN])
    # Invert to find survivors
    svr = .!ext
    pyplot(dpi=200)
    # Plot all the populations
    p1 = plot(yaxis=:log10,ylabel="Population (# cells)",ylims=(1e-5,Inf),xlabel="Time (s)")
    for i = 1:totN
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0)
        plot!(p1,T[inds],C[inds,i],label="")
    end
    savefig(p1,"Output/Fig3/all_pops.png")
    return(nothing)
end

@time figure3(111,500,1)
