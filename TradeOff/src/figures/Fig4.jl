# Script to plot figure 4
using TradeOff
using Plots
using JLD
using ColorSchemes
using LaTeXStrings
using Plots.PlotMeasures
import PyPlot

function sub_prob(M::Int64,R::Int64,subs::Float64)
    # Check that final waste product hasn't been generated
    if subs >= M - 1
        no_sub = M - 1
    # Also check that it hasn't somehow gone negative
    elseif subs <= 0.0
        no_sub = 0.0
    else
        no_sub = subs
    end
    # Calculate probability
    P = 1 - (1-(no_sub)/(M-1))^R
    return(P)
end

function figure4(rps::Int64,ims::Int64,sim_type::Int64,sim_type2::Int64)
    println("Compiled")
    # Extract other simulation parameters from the function
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Read in appropriate files
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    # Find file name to load in
    sfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/RunStats$(ims)Ims.jld"
    # Check it actually exists
    if ~isfile(sfile)
        error("missing stats file for $(ims) immigrations simulations")
    end
    # Read in relevant data
    ps = load(pfile,"ps")
    # Now load out the times, and number of trajectories
    times = load(sfile,"times")
    no_sims = load(sfile,"no_sims")
    no_via = load(sfile,"no_via")
    no_rs = load(sfile,"no_rs")
    # Load in averages
    mn_sbs = load(sfile,"mn_sbs")
    mn_via_R = load(sfile,"mn_via_R")
    mn_ηs_R = load(sfile,"mn_ηs_R")
    mn_KS_R = load(sfile,"mn_KS_R")
    # Load in standard deviations
    sd_sbs = load(sfile,"sd_sbs")
    sd_via_R = load(sfile,"sd_via_R")
    sd_ηs_R = load(sfile,"sd_ηs_R")
    sd_KS_R = load(sfile,"sd_KS_R")
    # Preallocate standard errors
    se_via_R = zeros(size(sd_via_R))
    se_ηs_R = zeros(size(sd_ηs_R))
    se_KS_R = zeros(size(sd_KS_R))
    # Calculate standard errors from this
    se_sbs = sd_sbs./sqrt.(no_sims)
    # Calculation (slightly) different in the viable case
    for i = 1:size(sd_via_R,1)
        se_via_R[i,:] = sd_via_R[i,:]./sqrt.(no_via)
        se_ηs_R[i,:] = sd_ηs_R[i,:]./sqrt.(no_rs[i,:])
        se_KS_R[i,:] = sd_KS_R[i,:]./sqrt.(no_rs[i,:])
    end
    # Preallocate probabilites
    mn_Ps = zeros(size(mn_via_R))
    up_Ps = zeros(size(mn_via_R))
    dw_Ps = zeros(size(mn_via_R))
    # Loop over, calculating the probability at each step
    for i = 1:size(mn_Ps,1)
        for j = 1:size(mn_Ps,2)
            # Calculate mean probability
            mn_Ps[i,j] = sub_prob(M,i,mn_sbs[j])
            # And also upper and lower bound
            up_Ps[i,j] = sub_prob(M,i,mn_sbs[j]+se_sbs[j]) .- mn_Ps[i,j]
            dw_Ps[i,j] = mn_Ps[i,j] .- sub_prob(M,i,mn_sbs[j]-se_sbs[j])
        end
    end
    println("Data read in")
    # Check if directory exists and if not make it
    if ~isdir("Output/Fig4")
        mkdir("Output/Fig4")
    end
    # Setup plotting
    pyplot(dpi=200)
    # Load in color scheme
    a = ColorSchemes.sunset.colors
    # Plot basic trade-off first
    p1 = plot(xlabel="Time (s)",ylabel="Number of strains",xlim=(-Inf,5e7))
    plot!(p1,title="High substrate free-energy case",legend=:bottomright)
    plot!(p1,times,mn_via_R[1,:],ribbon=se_via_R[1,:],label="R=1",color=a[1])
    plot!(p1,times,mn_via_R[3,:],ribbon=se_via_R[3,:],label="R=3",color=a[2])
    plot!(p1,times,mn_via_R[5,:],ribbon=se_via_R[5,:],label="R=5",color=a[3])
    plot!(p1,times,mn_via_R[7,:],ribbon=se_via_R[7,:],label="R=7",color=a[4])
    # Add annotation
    px, py = annpos([0.0;5e7],[0.0;7.6],0.05,0.0)
    annotate!(p1,px,py,text("A",17,:black))
    savefig(p1,"Output/Fig4/AvViaReacsTime.png")
    # Now do probability plot
    p2 = plot(xlabel="Time (s)",ylabel="Probability of no usable substrate",xlim=(-Inf,5e7))
    plot!(p2,title="Chance of species finding no usable substrates",legend=false)
    plot!(p2,times,1 .-mn_Ps[1,:],ribbon=(up_Ps[1,:],dw_Ps[1,:]),label="R=1",color=a[1])
    plot!(p2,times,1 .-mn_Ps[3,:],ribbon=(up_Ps[3,:],dw_Ps[3,:]),label="R=3",color=a[2])
    plot!(p2,times,1 .-mn_Ps[5,:],ribbon=(up_Ps[5,:],dw_Ps[5,:]),label="R=5",color=a[3])
    plot!(p2,times,1 .-mn_Ps[7,:],ribbon=(up_Ps[7,:],dw_Ps[7,:]),label="R=7",color=a[4])
    # Add annotation
    px, py = annpos([0.0;5e7],[0.0;1.1],0.05,0.0)
    annotate!(p2,px,py,text("B",17,:black))
    savefig(p2,"Output/Fig4/ProbSubTime.png")
    # Plot trade-off for η
    p3 = plot(xlabel="Time (s)",ylabel="Average eta value",xlim=(-Inf,5e7))
    plot!(p3,title="Variation of key parameters",legend=false)
    plot!(p3,times,mn_ηs_R[1,:],ribbon=se_ηs_R[1,:],label="R=1",color=a[1])
    plot!(p3,times,mn_ηs_R[3,:],ribbon=se_ηs_R[3,:],label="R=3",color=a[2])
    plot!(p3,times,mn_ηs_R[5,:],ribbon=se_ηs_R[5,:],label="R=5",color=a[3])
    plot!(p3,times,mn_ηs_R[7,:],ribbon=se_ηs_R[7,:],label="R=7",color=a[4])
    # Add annotation
    px, py = annpos([0.0;5e7],[1.2;3.95],0.05,0.0)
    annotate!(p3,px,py,text("C",17,:black))
    # Define box for inset here
    box = (1,bbox(0.4,0.15,0.4,0.3,:bottom,:left))
    Ks = L"K_S"
    e7 = L"10^7"
    em3 = L"10^{-3}"
    # Plot other trade-off into the inset
    plot!(p3,times/1e7,mn_KS_R[1,:]*1000.0,ribbon=se_KS_R[1,:]*1000.0,color=a[1],label="",inset_subplots=box,subplot=2)
    plot!(p3,times/1e7,mn_KS_R[3,:]*1000.0,ribbon=se_KS_R[3,:]*1000.0,color=a[2],label="",subplot=2)
    plot!(p3,times/1e7,mn_KS_R[5,:]*1000.0,ribbon=se_KS_R[5,:]*1000.0,color=a[3],label="",subplot=2)
    plot!(p3,times/1e7,mn_KS_R[7,:]*1000.0,ribbon=se_KS_R[7,:]*1000.0,color=a[4],label="",subplot=2)
    plot!(p3,xlabel="Time ($(e7) s)",ylabel="$(Ks) ($em3)",xlim=(-Inf,5.0),grid=false,subplot=2)
    savefig(p3,"Output/Fig4/AvEtaperReacTime.png")
    # Extract simulation parameters for the other case from the function
    Np, Nt, M, d, μrange = sim_paras(sim_type2)
    # Find file name to load in
    sfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/RunStats$(ims)Ims.jld"
    # Check it actually exists
    if ~isfile(sfile)
        error("missing stats file for $(ims) immigrations simulations")
    end
    # Now load out the times, and number of trajectories
    times = load(sfile,"times")
    no_via = load(sfile,"no_via")
    # Load in averages
    mn_via_R = load(sfile,"mn_via_R")
    # Load in standard deviations
    sd_via_R = load(sfile,"sd_via_R")
    # Preallocate standard errors
    se_via_R = zeros(size(sd_via_R))
    # Calculation (slightly) different in the viable case
    for i = 1:size(sd_via_R,1)
        se_via_R[i,:] = sd_via_R[i,:]./sqrt.(no_via)
    end
    p4 = plot(xlabel="Time (s)",ylabel="Number of strains",xlim=(-Inf,5e7))
    plot!(p4,title="Low substrate free-energy case",legend=false)
    plot!(p4,times,mn_via_R[1,:],ribbon=se_via_R[1,:],label="R=1",color=a[1])
    plot!(p4,times,mn_via_R[3,:],ribbon=se_via_R[3,:],label="R=3",color=a[2])
    plot!(p4,times,mn_via_R[5,:],ribbon=se_via_R[5,:],label="R=5",color=a[3])
    plot!(p4,times,mn_via_R[7,:],ribbon=se_via_R[7,:],label="R=7",color=a[4])
    # Add annotation
    px, py = annpos([0.0;5e7],[0.0;9.25],0.05,0.0)
    annotate!(p4,px,py,text("D",17,:black))
    savefig(p4,"Output/Fig4/LowFreeEnergy.png")
    # Plot all graphs as a single figure
    pt = plot(p1,p3,p2,p4,layout=4,size=(1200,800),margin=5.0mm)
    savefig(pt,"Output/Fig4/figure4.png")
    return(nothing)
end

@time figure4(250,500,1,2)
