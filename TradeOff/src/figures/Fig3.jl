# Script to plot figure 3
using TradeOff
using Plots
using JLD
using ColorSchemes
using LaTeXStrings
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
    lfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/SurvTimes$(ims)Ims.jld"
    if ~isfile(lfile)
        error("$(ims) immigrations simulation is missing a long term survivors file")
    end
    # Read in relevant data
    sTs = load(lfile,"sTs")
    # Now look at snapshot data
    sfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/SnapDataStats$(ims)Ims.jld"
    if ~isfile(sfile)
        error("$(ims) immigrations simulation is missing an output file")
    end
    # Save growth probabilities
    gp = load(sfile,"gp")
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
    # Load in final ϕR values
    all_fin_ϕRs = load(tfile,"all_fin_ϕRs")
    # Load in averages
    mn_svt = load(tfile,"mn_svt")
    mn_pop = load(tfile,"mn_pop")
    mn_shD = load(tfile,"mn_shD")
    # Load in standard deviations
    sd_svt = load(tfile,"sd_svt")
    sd_pop = load(tfile,"sd_pop")
    sd_shD = load(tfile,"sd_shD")
    # Calculate standard errors
    se_svt = sd_svt./sqrt.(no_sims)
    se_pop = sd_pop./sqrt.(no_sims)
    se_shD = sd_shD./sqrt.(no_sims)
    println("Data read in")
    # Check if directory exists and if not make it
    if ~isdir("Output/Fig3")
        mkdir("Output/Fig3")
    end
    # Setup plotting
    pyplot(dpi=200)
    # Define latex commands
    e7 = L"10^7"
    e13 = L"10^{13}"
    # Load in color scheme
    a = ColorSchemes.tab10.colors
    # Make first plot
    p1 = plot(ylabel="Total population ($(e13) cells)",xlim=(-Inf,5e7),legend=:topleft)
    plot!(p1,t_times,mn_pop/1e13,ribbon=se_pop/1e13,label="Population",color=a[1],ylim=(-Inf,4.0))
    # Define box for inset here
    box = (1,bbox(0.4,0.35,0.4,0.3,:bottom,:left))
    # Add histogram in as an insert
    histogram!(p1,all_fin_ϕRs,nbins=100,color=a[2],label="",inset_subplots=box,subplot=2)
    plot!(p1,xlabel="Final ribosome fraction ($(L"\phi_R"))",ylabel="Number of survivors",subplot=2)
    plot!(p1,guidefontsize=9,grid=false,subplot=2)
    # Twin the x-axis
    pt = twinx(p1)
    # Ensure same limits are used
    plot!(pt,xlim=(-Inf,5e7),ylim=(-Inf,2.9),ylabel="Shannon diversity")
    # Then plot the Shannon diversity
    p1 = plot!(pt,t_times,mn_shD,ribbon=se_shD,label="Diversity",color=a[3],legend=:topright)
    # Add annotation
    px, py = annpos([0.0;5e7],[0.0;3.9],0.05,0.0)
    annotate!(p1,px,py,text("A",17,:black))
    savefig(p1,"Output/Fig3/SumStats.png")
    # Now make the second plot
    p2 = plot(xlabel="Time (s)",ylabel="Number of species",xlim=(-Inf,5e7),legend=:topleft)
    plot!(p2,t_times,mn_svt,ribbon=se_svt,label="Species",color=a[4])
    # Add annotation
    px, py = annpos([0.0;5e7],[0.0;35.0],0.05,0.0)
    annotate!(p2,px,py,text("B",17,:black))
    # Define box for inset here
    box = (1,bbox(0.4,0.4,0.4,0.3,:bottom,:left))
    # Add histogram in as an insert
    histogram!(p2,sTs/1e7,nbins=100,color=a[5],label="",inset_subplots=box,subplot=2)
    plot!(p2,xlabel="Time of immigration ($(e7) s)",ylabel="Number of survivors",subplot=2)
    plot!(p2,guidefontsize=9,grid=false,subplot=2)
    # Twin the x-axis
    pt = twinx(p2)
    # Ensure same limits are used
    plot!(pt,xlim=(-Inf,5e7),ylabel="Proportion of immigrants that can grow")
    # Then plot the proportion feasible
    p2 = plot!(pt,s_times,gp,label="Invasibility",color=a[6],legend=:topright)
    savefig(p2,"Output/Fig3/Invasibility.png")
    # Now want to make a plot incorperating both previous plots
    pt = plot(p1,p2,layout=(2,1),size=(750,800),margin=5.0mm,rightmargin=20mm)
    savefig(pt,"Output/Fig3/figure3.png")
    return(nothing)
end

@time figure3(250,500,1)
