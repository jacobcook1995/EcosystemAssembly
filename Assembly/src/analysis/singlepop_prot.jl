# A script to analyse the proteome model for a single population.
using Assembly
using Plots
using LaTeXStrings
import PyPlot

# Just want to set up and run a single population
function singpop()
    println("Successfully compiled.")
    # Simple test data set
    d = 5e-3 # death rate
    κ = [100.0,0.0] # Metabolite supply rate
    δ = 1.0*ones(2) # Metabolite dilution rate
    KS = 0.1 # Saturation constant
    kr = 10.0 # Reversibility factor
    k = 1.0 # matches qm from previously
    E = 1.0 # FIXED FOR NOW BUT SHOULD BE ADJUSTED LATER
    ai = 5.0 # initial energy level
    Ni = 100.0 # initial population
    ρ = 1e-7
    # Initialise parameter set
    ps = initialise_prot()
    # Choose initial protein fractions
    ϕ = [0.275,0.275,0.45] # Again this should shift
    pa = make_var_prot(ps,ϕ)
    # Choose simulation time
    Tmax = 2500.0
    # Then run simulation
    C, T = prot_simulate(ps,Tmax,ai,Ni,pa,d,κ,δ,KS,kr,k,ρ)

    # Do plotting
    pyplot(dpi=200)
    plot(T,C[:,1],xlabel="Time",label="",ylabel="Population")
    savefig("Output/testPop.png")
    plot(T,C[:,2],xlabel="Time",label="",ylabel="Cell energy conc")
    savefig("Output/testEng.png")
    plot(T,C[:,3:4],xlabel="Time",label=["Substrate" "Waste"],ylabel="Concentration")
    savefig("Output/testCon.png")
    return(nothing)
end

@time singpop()
