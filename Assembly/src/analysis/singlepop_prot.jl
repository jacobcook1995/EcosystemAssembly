# A script to analyse the proteome model for a single population.
using Assembly
using Plots
using LaTeXStrings
import PyPlot

# Just want to set up and run a single population
function singpop()
    println("Successfully compiled.")
    # Simple test data set
    d = 5e-3 # death rate # THIS ONE IS PRETTY ARBITARY, ADJUST TO FIT. => MAYBE LINK TO δ???
    ai = 5.0 # initial energy level
    Ni = 100.0 # initial population
    ρ = 1.0
    # Initialise parameter set
    ps = initialise_prot()
    # Choose initial protein fractions
    ϕ = [0.275,0.275,0.45] # Again this should shift
    pa = make_var_prot(ps,ϕ)
    # Choose simulation time
    Tmax = 10000.0
    # Then run simulation
    C, T = prot_simulate(ps,Tmax,ai,Ni,pa,d,ρ)

    λa = zeros(length(T))
    for i = 1:length(T)
        λa[i] = λs(C[i,2],ϕ[1],ps)
    end
    # Do plotting
    pyplot(dpi=200)
    plot(T,C[:,1],xlabel="Time",label="",ylabel="Population")
    savefig("Output/testPop.png")
    plot(T,C[:,2],xlabel="Time",label="",ylabel="Cell energy conc")
    savefig("Output/testEng.png")
    plot(T,C[:,3:4],xlabel="Time",label=["Substrate" "Waste"],ylabel="Concentration")
    savefig("Output/testCon.png")
    plot(T,λa)
    savefig("Output/test.png")
    return(nothing)
end

@time singpop()
