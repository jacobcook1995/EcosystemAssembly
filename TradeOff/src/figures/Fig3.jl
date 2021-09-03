# Script to plot figure 3
using TradeOff
using Plots
using JLD
using Plots.PlotMeasures
import PyPlot

function figure3(rN::Int64,ims::Int64,sim_type::Int64)
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
    # Plot populations that survive to the end
    p2 = plot(yaxis=:log10,ylabel="Population (# cells)",ylims=(1e-5,Inf),xlims=(0.0,Inf),xlabel="Time (s)")
    for i = 1:totN
        if svr[i] == true
            # Find and eliminate zeros so that they can be plotted on a log plot
            inds = (C[:,i] .> 0)
            plot!(p2,T[inds],C[inds,i],label="")
        end
    end
    # Add annotation
    px, py = annpos([0.0;5e7],[1e-5;5e12],0.05,0.0)
    annotate!(p2,px,py,text("A",17,:black))
    savefig(p2,"Output/Fig3/surv_pops.png")
    # Plot ribosome fractions of populations that survive to the end
    p3 = plot(ylabel="Ribosome fraction",xlims=(0.0,Inf),xlabel="Time (s)")
    for i = 1:totN
        if svr[i] == true
            plot!(p3,T,C[:,2*totN+ps.M+i],label="")
        end
    end
    # Add annotation
    px, py = annpos([0.0;5e7],[0.02;0.215],0.05,0.0)
    annotate!(p3,px,py,text("B",17,:black))
    savefig(p3,"Output/Fig3/surv_fracs.png")
    # Now want to make a plot incorperating both previous plots
    pt = plot(p2,p3,layout=(2,1),size=(700,800),margin=5.0mm)
    savefig(pt,"Output/Fig3/figure3.png")
    return(nothing)
end

@time figure3(111,500,1)
