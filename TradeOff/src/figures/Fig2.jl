# Script to plot elements needed for figure 1
using TradeOff
using Plots
using JLD
import PyPlot

# Function to plot the efficency with changing ribsome fraction
function early_immigration_dyns(rN::Int64,ims::Int64,sim_type::Int64)
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
    if ~isdir("Output/Fig2")
        mkdir("Output/Fig2")
    end
    # Find total number of strains
    totN = length(micd)
    # Set maximum time to plot to
    Tmax = 5e5
    pyplot(dpi=200)
    # Plot all the populations
    p1 = plot(yaxis=:log10,ylabel="Population (# cells)",ylims=(1e-5,Inf),xlabel="Time (s)")
    for i = 1:totN
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0) .& (T .<= Tmax)
        plot!(p1,T[inds],C[inds,i],label="")
    end
    savefig(p1,"Output/Fig2/all_pops.png")
    return(nothing)
end

@time early_immigration_dyns(3,500,1)
