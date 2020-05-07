# Script that reads in data generated using repeat and then makes plots using it
using Assembly
using Plots
using JLD
import PyPlot

# Plotting function
function plot_repeat()
    if length(ARGS) == 0
        error("NEED TO PROVIDE NAME FOR INPUT DATA.")
    end
    # First read in parameter set
    ps = load("Paras/ps$(ARGS[1]).jld","ps")
    # Make reaction matrix
    rM = zeros(Int64,ps.M,ps.M)
    for i = 1:length(ps.reacs)
        rM[ps.reacs[i].Prd,ps.reacs[i].Rct] += 1
    end
    pyplot(dpi=200)
    heatmap(1:size(rM,1),1:size(rM,2),rM,colorbar=:none)
    plot!(xlabel="Reactant",ylabel="Product")
    plot!(0:100,0:100,color=:white,label="")
    savefig("Output/RandReacs$(ARGS[1]).png")
    # Now read in dynamics data
    C = load("Output/Dyn$(ARGS[1]).jld","C")
    T = load("Output/Dyn$(ARGS[1]).jld","T")
    t = load("Output/Dyn$(ARGS[1]).jld","t")
    # Need to recalculate N
    N = ps.N + length(t) - 1
    # Now plot the dynamics
    plot(T,C[:,1:N],label="",xlabel="Time",ylabel="Population")
    vline!(t,color=:red,label="",style=:dot)
    savefig("Output/Pops$(ARGS[1]).png")
    plot(T,C[:,N+1:end],label="")
    vline!(t,color=:red,label="",style=:dot)
    savefig("Output/Concs$(ARGS[1]).png")
    return(nothing)
end

@time plot_repeat()
