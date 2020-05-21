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
    # Choose initial protein fractions
    ϕ = [0.275,0.275,0.45] # Again this should shift
    pa = make_var_prot(ps,ϕ)
    # Choose simulation time
    Tmax = 1000000.0
    # Then run simulation
    C, T = prot_simulate(ps,Tmax,ai,Ni,pa)

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
    s1 = L"s^{-1}"
    plot(T,λa,xlabel="Time",label="",ylabel="Growth rate $(s1)")
    savefig("Output/testGrowth.png")
    return(nothing)
end

# Similar function as above except that it tries to find optimal values for ϕ_M
function singpop_opt()
    println("Successfully compiled.")
    # Simple test data set
    ai = 5.0 # initial energy level
    Ni = 100.0 # initial population
    # Initialise parameter set
    ps = initialise_prot(true)
    # Pick an initial condition to optimise for (assume fixed environment)
    S = 100.0
    P = [collect(0.0:1.0:99.0); collect(99.0:0.1:99.9); collect(99.90:0.01:99.99); 99.99999999]
    θs = zeros(length(P))
    # Make vector of theta values
    for i = 1:length(P)
        θs[i] = θ(S,P[i],ps.T,ps.η,ps.r.ΔG0)
    end
    # Now find optimal ribosome fractions and growth rates for each one
    ϕR = zeros(length(P))
    λs = zeros(length(P))
    ao = zeros(length(P))
    # Housekeeping fraction is fixed throughout
    ϕH = 0.45
    for i = 1:length(P)
        # Use function to find optimal ribosome fraction
        ϕ, λs[i], ao[i] = optimise_ϕ(S,P[i],ps,ϕH)
        ϕR[i] = ϕ[1]
    end
    pyplot(dpi=200)
    plot(θs,ϕR,label="",xlabel=L"θ",ylabel=L"ϕ_R")
    savefig("Output/RibosomeFrac.png")
    plot(θs,λs,label="",xlabel=L"θ",ylabel="Optimal growth rate")
    savefig("Output/GrowthRate.png")
    plot(θs,ao,label="",xlabel=L"θ",ylabel="ATP per cell")
    savefig("Output/EnergyConc.png")
    return(nothing)
end

# Run a single population multiple times and plot a scatter graph of the results
function singpop_scat()
    println("Successfully compiled.")
    # Simple test data set
    ai = 5.0 # initial energy level
    Ni = 100.0 # initial population
    # Initialise parameter set
    ps = initialise_prot(false)
    # Choose simulation time
    Tmax = 500000.0
    # Then run multiple simulations
    a, J = prot_simulate_mult(ps,ai,Ni,Tmax)
    pyplot(dpi=200)
    m = L"^{-1}"
    plot(xlabel="Energy acquisition rate ATP cell$(m) s$m",ylabel="ATP per cell")
    for i = 1:9
        # Tom used a log plot but I think this will obsurce too much at the moment
        scatter!(J[i,:],a[i,:],label="$i")
    end
    savefig("Output/test.png")
    return(nothing)
end

@time singpop()
