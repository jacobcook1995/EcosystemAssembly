# A script to analyse single population growth to test and parameterise the new model
using TradeOff
using Random
using Plots
import PyPlot

# Simply want to test if the new function I've written works
function test()
    println("Compiled and starting simulation")
    # Two metabolites 1 reaction for the sake of simplicity
    M = 2
    O = 1
    # First generate parameter set
    ps = initialise(M,O)
    # Only one reaction so upper and lower bound have to be one
    Rl = 1
    Ru = 1
    # Choose new parameter values
    μ = 1e-1
    ω = 0.5
    KΩ = 1e9
    Kχ = 1e9
    # Then make a microbes
    mic = new_mic(M,μ,ω,KΩ,Kχ,Rl,Ru)
    # Set reasonable time window
    Tmax = 1e6
    # Choose initial condition
    pop = 1000.0
    conc = 1e-5
    as = 1e5
    ϕs = 0.128
    # Then simulate
    C, T = sing_pop(ps,pop,conc,as,ϕs,mic,Tmax)
    # Finally plot this
    pyplot(dpi=200)
    # Only 1 strain at the moment
    totN = 1
    # Plot all the populations
    p1 = plot(yaxis=:log10,ylabel="Population (# cells)",ylims=(1e-5,Inf))
    for i = 1:totN
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0)
        plot!(p1,T[inds],C[inds,i],label="")
    end
    savefig(p1,"Output/pops.png")
    # Plot all the concentrations
    p2 = plot(yaxis=:log10,ylabel="Concentration")#,ylims=(1e-15,Inf))
    for i = 1:ps.M
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,totN+i] .> 0)
        plot!(p2,T[inds],C[inds,totN+i],label="")
    end
    savefig(p2,"Output/concs.png")
    plot(T,C[:,(totN+ps.M+1):(2*totN+ps.M)],label="")
    savefig("Output/as.png")
    plot(T,C[:,(2*totN+ps.M+1):end],label="")
    savefig("Output/fracs.png")
    return(nothing)
end

@time test()
