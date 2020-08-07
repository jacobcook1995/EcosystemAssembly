# A script to form stable communities using the full model that then saves them for later use
using Assembly
using Plots
using LaTeXStrings
import PyPlot

# function to test that the new stuff I'm writing actually works
function test()
    println("Compiled!")
    # Assume that half saturation occurs at a quarter κ/δ
    KS = (1/4)*5.5e-3
    # From wikipedia an average enzyme has a k to KS ratio of 10^5 M^-1 s^-1
    # This would give us a k of 137.5, sensible to assume an above average rate
    # Though should be reduced by the fact we include uptake as well as metabolism
    # Choosing k = 500 means we match the maximum glucose uptake rate seen in Natarajan et al (2000)
    # of 3*10^7 molecules per second.
    # The above is a sensible argument but 1.0 gives a more reasonable ATP concentration.
    kc = 1.0
    # The reversibility factor remains the same as previously
    kr = 10.0
    # Assume microbes have 2 reactions each
    mR = 2.0
    sdR = 0.0
    # Case of 8 metabolites and 1 strain
    N = 1
    M = 8
    # Use formula to find how many reactions this implies
    O = 2*M - 3
    good = false
    # Make parameter set
    ps = 0
    while good == false
        ps = initialise(N,M,O,mR,sdR,kc,KS,kr)
        # This step is to ensure that the first metabolite can be broken down
        if any((ps.reacs[ps.mics[1].Reacs].↦:Rct) .== 1)
            good = true
        end
    end
    # Set time long enough to see dynamics
    Tmax = 1000000.0
    # Fairly arbitary inital conditions
    pop = ones(N)
    conc = zeros(M)
    as = 1e5*ones(N)
    ϕs = 0.2*ones(N)
    C, T = full_simulate(ps,Tmax,pop,conc,as,ϕs)
    pyplot(dpi=200)
    # Find and eliminate zeros so that they can be plotted on a log plot
    x = zeros(length(T),N)
    for i = 1:N
        x[:,i] = T
    end
    y = C[:,1:N]
    inds = (y .> 0)
    plot(x[inds],y[inds],xlabel="Time",label="",ylabel="Population",yaxis=:log10)
    savefig("Output/PopvsTime.png")
    plot(T,C[:,N+1:N+M],xlabel="Time",label="",ylabel="Concentration")
    savefig("Output/MetabolitevsTime.png")
    plot(T,C[:,N+M+1:2*N+M],xlabel="Time",label="",ylabel="Cell energy conc")
    savefig("Output/EnergyvsTime.png")
    plot(T,C[:,2*N+M+1:end],xlabel="Time",label="",ylabel=L"\phi_R")
    savefig("Output/FractionvsTime.png")
    return(nothing)
end

@time test()
