# Script to plot figure 5
using TradeOff
using Plots
using JLD
using ColorSchemes
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

function figure5(rps::Int64,ims::Int64,sim_type::Int64)
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
    # Check if directory exists and if not make it
    if ~isdir("Output/Fig5")
        mkdir("Output/Fig5")
    end
    # Setup plotting
    pyplot(dpi=200)
    # Load in color scheme
    a = ColorSchemes.sunset.colors
    # Plot basic trade-off first
    p1 = plot(xlabel="Time (s)",ylabel="Number of strains",xlim=(-Inf,5e7),legend=:bottomright)
    plot!(p1,times,mn_via_R[1,:],ribbon=se_via_R[1,:],label="R=1",color=a[1])
    plot!(p1,times,mn_via_R[3,:],ribbon=se_via_R[3,:],label="R=3",color=a[2])
    plot!(p1,times,mn_via_R[5,:],ribbon=se_via_R[5,:],label="R=5",color=a[3])
    plot!(p1,times,mn_via_R[7,:],ribbon=se_via_R[7,:],label="R=7",color=a[4])
    # Add annotation
    px, py = annpos([0.0;5e7],[0.0;7.0],0.05,0.0)
    annotate!(p1,px,py,text("A",17,:black))
    savefig(p1,"Output/Fig5/AvViaReacsTime.png")
    # Now do probability plot
    p2 = plot(xlabel="Time (s)",ylabel="Probability of usable substrate",xlim=(-Inf,5e7),legend=false)
    plot!(p2,times,mn_Ps[1,:],ribbon=(dw_Ps[1,:],up_Ps[1,:]),label="R=1",color=a[1])
    plot!(p2,times,mn_Ps[3,:],ribbon=(dw_Ps[3,:],up_Ps[3,:]),label="R=3",color=a[2])
    plot!(p2,times,mn_Ps[5,:],ribbon=(dw_Ps[5,:],up_Ps[5,:]),label="R=5",color=a[3])
    plot!(p2,times,mn_Ps[7,:],ribbon=(dw_Ps[7,:],up_Ps[7,:]),label="R=7",color=a[4])
    # Add annotation
    px, py = annpos([0.0;5e7],[0.0;1.1],0.05,0.0)
    annotate!(p2,px,py,text("B",17,:black))
    savefig(p2,"Output/Fig5/ProbSubTime.png")
    # Plot trade-off for η
    p3 = plot(xlabel="Time (s)",ylabel="Average eta value",xlim=(-Inf,5e7),legend=false)
    plot!(p3,times,mn_ηs_R[1,:],ribbon=se_ηs_R[1,:],label="R=1",color=a[1])
    plot!(p3,times,mn_ηs_R[3,:],ribbon=se_ηs_R[3,:],label="R=3",color=a[2])
    plot!(p3,times,mn_ηs_R[5,:],ribbon=se_ηs_R[5,:],label="R=5",color=a[3])
    plot!(p3,times,mn_ηs_R[7,:],ribbon=se_ηs_R[7,:],label="R=7",color=a[4])
    # Add annotation
    px, py = annpos([0.0;5e7],[1.2;3.75],0.05,0.0)
    annotate!(p3,px,py,text("C",17,:black))
    savefig(p3,"Output/Fig5/AvEtaperReacTime.png")
    # Plot trade-off for the saturation constant
    p4 = plot(xlabel="Time (s)",ylabel="Average half saturation constant",xlim=(-Inf,5e7),legend=false)
    plot!(p4,times,mn_KS_R[1,:],ribbon=se_KS_R[1,:],label="R=1",color=a[1])
    plot!(p4,times,mn_KS_R[3,:],ribbon=se_KS_R[3,:],label="R=3",color=a[2])
    plot!(p4,times,mn_KS_R[5,:],ribbon=se_KS_R[5,:],label="R=5",color=a[3])
    plot!(p4,times,mn_KS_R[7,:],ribbon=se_KS_R[7,:],label="R=7",color=a[4])
    # Add annotation
    px, py = annpos([0.0;5e7],[0.0005;0.00215],0.05,0.0)
    annotate!(p4,px,py,text("D",17,:black))
    savefig(p4,"Output/Fig5/AvKSperReacTime.png")
    # Plot all graphs as a single figure
    pt = plot(p1,p3,p2,p4,layout=4,size=(1200,800),margin=5.0mm)
    savefig(pt,"Output/Fig5/figure5.png")
    return(nothing)
end

@time figure5(250,500,1)
