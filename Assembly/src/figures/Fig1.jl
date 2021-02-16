# Script to plot elements needed for figure 1
using Assembly
using Plots
using JLD
import PyPlot

# function to plot population dynamics
function popdyn(Rl::Int64,Ru::Int64,syn::Bool,Nr::Int64,Ns::Int64,en::String)
    println("Compiled!")
    # Read in relevant files
    pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/ParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(Nr)Ns$(Ns).jld"
    if ~isfile(pfile)
        error("run $(Nr) is missing a parameter file")
    end
    ofile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/OutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(Nr)Ns$(Ns).jld"
    if ~isfile(ofile)
        error("run $(Nr) is missing an output file")
    end
    efile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/ExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(Nr)Ns$(Ns).jld"
    if ~isfile(efile)
        error("run $(Nr) is missing an extinct file")
    end
    ps = load(pfile,"ps")
    C = load(ofile,"C")
    T = load(ofile,"T")
    out = load(ofile,"out")
    ded = load(efile,"ded")
    # Find maximum time
    Tmax = T[end]
    # Now move onto plotting
    pyplot()
    theme(:wong2,dpi=200)
    # Plot all the populations
    plot(title="Population dynamics",yaxis=:log10,xlabel="Time (s)",ylabel="Population (# cells)")
    for i = 1:Ns
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0) .& (T .< Tmax/4)
        plot!(T[inds],C[inds,i],label="")
    end
    savefig("Output/Fig1/fullpops.png")
    return(nothing)
end

@time popdyn(1,5,true,78,250,"l")
