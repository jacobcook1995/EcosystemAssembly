# A script to analyse the proteome model for a single population.
using Assembly
using Plots
using LaTeXStrings
import PyPlot

# Just want to set up and run a single population
function singpop()
    println("Successfully compiled.")
    # Simple test data set
    ai = 5.0 # initial energy level
    Ni = 100.0 # initial population
    # Initialise parameter set
    ps = initialise_prot(false)
    # Choose simulation time
    Tmax = 1000000.0
    # Then run simulation
    C, T = prot_simulate(ps,Tmax,ai,Ni)
    # Now calculate growth rates and proteome fractions
    λa = zeros(length(T))
    ϕR = zeros(length(T))
    for i = 1:length(T)
        ϕR[i] = ϕ_R(C[i,2],ps)
        λa[i] = λs(C[i,2],ϕR[i],ps)
    end
    # Do plotting
    pyplot(dpi=200)
    plot(T,log10.(C[:,1]),xlabel="Time",label="",ylabel="Population")
    savefig("Output/PopvsTime.png")
    plot(T,C[:,2],xlabel="Time",label="",ylabel="Cell energy conc")
    savefig("Output/EnergyvsTime.png")
    plot(T,C[:,3:4],xlabel="Time",label=["Substrate" "Waste"],ylabel="Concentration")
    savefig("Output/MetabolitevsTime.png")
    s1 = L"s^{-1}"
    plot(T,λa,xlabel="Time",label="",ylabel="Growth rate $(s1)")
    savefig("Output/GrowthvsTime.png")
    plot(T,ϕR,xlabel="Time",label="",ylabel=L"\phi_R")
    plot!(T,C[:,5],label="")
    savefig("Output/FractionvsTime.png")
    return(nothing)
end

@time singpop()
