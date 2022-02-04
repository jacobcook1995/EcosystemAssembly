# Script to plot figure 7 (which shows situation without immigration)
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

# Function to plot the 7th figure (e.g. development without immigration)
function figure7(rps::Int64)
    println("Compiled")
    # No immigration situation (aka sim type 5)
    ims = 0
    sim_type = 5
    tk = "NoImm"
    # Extract other simulation parameters from the function
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Read in appropriate files
    pfile = "Output/$(tk)$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    # Find file name to load in
    sfile = "Output/$(tk)$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/RunStats$(ims)Ims.jld"
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
    # Load in averages
    mn_sbs = load(sfile,"mn_sbs")
    mn_via_R = load(sfile,"mn_via_R")
    mn_via_ω = load(sfile,"mn_via_ω")
    mn_via_ϕR = load(sfile,"mn_via_ϕR")
    mn_via_η = load(sfile,"mn_via_η")
    # Load in standard deviations
    sd_sbs = load(sfile,"sd_sbs")
    sd_via_R = load(sfile,"sd_via_R")
    sd_via_ω = load(sfile,"sd_via_ω")
    sd_via_ϕR = load(sfile,"sd_via_ϕR")
    sd_via_η = load(sfile,"sd_via_η")
    # Preallocate standard errors
    se_via_R = zeros(size(sd_via_R))
    # Calculate standard errors from this
    se_sbs = sd_sbs./sqrt.(no_sims)
    se_via_ω = sd_via_ω./sqrt.(no_via)
    se_via_ϕR = sd_via_ϕR./sqrt.(no_via)
    se_via_η = sd_via_η./sqrt.(no_via)
    # Calculation (slightly) different in the viable case
    for i = 1:size(sd_via_R,1)
        se_via_R[i,:] = sd_via_R[i,:]./sqrt.(no_via)
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
    if ~isdir("Output/Fig7")
        mkdir("Output/Fig7")
    end
    # Setup plotting
    pyplot(dpi=200)
    # Load in color scheme
    a = ColorSchemes.sunset.colors
    # Plot basic trade-off first
    p1 = plot(xlabel="Time (s)",ylabel="Number of species",xlim=(-Inf,2.5e6))
    plot!(p1,title="Number of reactions with time",legend=:bottomright)
    plot!(p1,times,mn_via_R[1,:],ribbon=se_via_R[1,:],label="R=1",color=a[1])
    plot!(p1,times,mn_via_R[3,:],ribbon=se_via_R[3,:],label="R=3",color=a[2])
    plot!(p1,times,mn_via_R[5,:],ribbon=se_via_R[5,:],label="R=5",color=a[3])
    plot!(p1,times,mn_via_R[7,:],ribbon=se_via_R[7,:],label="R=7",color=a[4])
    # Add annotation
    px, py = annpos([0.0;2.5e6],[0.0;10.0],0.05,0.0)
    annotate!(p1,px,py,text("A",17,:black))
    savefig(p1,"Output/Fig7/AvViaReacsTime.png")
    # Now do probability plot
    p2 = plot(xlabel="Time (s)",ylabel="Probability of no usable substrate",xlim=(-Inf,2.5e6))
    plot!(p2,title="Chance of species finding no usable substrates",legend=false)
    plot!(p2,times,1 .-mn_Ps[1,:],ribbon=(up_Ps[1,:],dw_Ps[1,:]),label="R=1",color=a[1])
    plot!(p2,times,1 .-mn_Ps[3,:],ribbon=(up_Ps[3,:],dw_Ps[3,:]),label="R=3",color=a[2])
    plot!(p2,times,1 .-mn_Ps[5,:],ribbon=(up_Ps[5,:],dw_Ps[5,:]),label="R=5",color=a[3])
    plot!(p2,times,1 .-mn_Ps[7,:],ribbon=(up_Ps[7,:],dw_Ps[7,:]),label="R=7",color=a[4])
    # Add annotation
    px, py = annpos([0.0;2.5e6],[0.0;1.1],0.05,0.0)
    annotate!(p2,px,py,text("B",17,:black))
    savefig(p2,"Output/Fig7/ProbSubTime.png")
    # Load in color scheme
    a = ColorSchemes.Dark2_4.colors
    # Make box for the subplots
    box = (1,bbox(0.5,0.65,0.4,0.3,:bottom,:left))
    # Define reused latex strings
    e6 = L"10^6"
    # Make plot objects
    p3 = plot(xlabel="Time (s)",xlim=(-Inf,2.5e6),ylim=(0.5,0.65),legend=:topleft)
    plot!(p3,title="Ribosome dynamics",ylabel="Maximum ribsome fraction factor ($(L"\omega"))")
    # Add inset box to plot other trade-off into
    plot!(p3,xlabel="Time ($(e6) s)",ylabel=L"\phi_R",inset_subplots=box,subplot=2,grid=false,legend=false)
    plot!(p3,subplot=2,xlim=(-Inf,2.5),ylim=(0.05,0.3),legend=false)
    # Then actually plot the relevant data
    plot!(p3,times,mn_via_ω,ribbon=se_via_ω,label="",color=a[1])
    plot!(p3,times/1e6,mn_via_ϕR,ribbon=se_via_ϕR,color=a[1],subplot=2)
    # Add annotation
    px, py = annpos([0.0;2.5e6],[0.0;0.654],0.1,0.0)
    annotate!(p3,px,py,text("C",17,:black))
    savefig(p3,"Output/Fig7/omega_with_t.png")
    # Make plot objects
    p4 = plot(xlabel="Times (s)",xlim=(-Inf,2.5e6),legend=:right,ylabel=L"\eta",title="ATP yield with time")
    # Plot the data to the relevant plot objects
    plot!(p4,times,mn_via_η,ribbon=se_via_η,color=a[1],label="")
    # Add annotation
    px, py = annpos([0.0;2.5e6],[0.0;3.385],0.05,0.0)
    annotate!(p4,px,py,text("D",17,:black))
    # Save figures to this directory
    savefig(p4,"Output/Fig7/Eta.png")
    # Plot all graphs as a single figure
    pt = plot(p1,p3,p2,p4,layout=4,size=(1200,800),margin=5.0mm)
    savefig(pt,"Output/Fig7/figure7.png")
    return(nothing)
end

@time figure7(250)
