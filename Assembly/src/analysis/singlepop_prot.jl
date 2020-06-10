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
    plot(T,C[:,1],xlabel="Time",label="",ylabel="Population")
    savefig("Output/PopvsTime.png")
    plot(T,C[:,2],xlabel="Time",label="",ylabel="Cell energy conc")
    savefig("Output/EnergyvsTime.png")
    plot(T,C[:,3:4],xlabel="Time",label=["Substrate" "Waste"],ylabel="Concentration")
    savefig("Output/MetabolitevsTime.png")
    s1 = L"s^{-1}"
    plot(T,λa,xlabel="Time",label="",ylabel="Growth rate $(s1)")
    savefig("Output/GrowthvsTime.png")
    plot(T,ϕR,xlabel="Time",label="",ylabel=L"\phi_R")
    savefig("Output/FractionvsTime.png")
    return(nothing)
end

# Similar function as above except that it tries to find optimal values for ϕ_M
function singpop_opt()
    println("Successfully compiled.")
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
    ϕR2 = zeros(length(P))
    λs = zeros(length(P))
    ao = zeros(length(P))
    for i = 1:length(P)
        # Use function to find optimal ribosome fraction
        ϕ, λs[i], ao[i] = optimise_ϕ(S,P[i],ps)
        ϕR[i] = ϕ[1]
        # Find what fraction I am using
        ϕR2[i] = ϕ_R(ao[i],ps)
    end
    pyplot(dpi=200)
    plot(θs,ϕR,label="",xlabel=L"θ",ylabel=L"ϕ_R")
    plot!(θs,ϕR2,label="")
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
    # Make lists
    ls = Array{Any,1}(undef,7)
    # 13 decline
    ls[1] = collect(4:12)
    # 14 decline
    ls[2] = collect(4:13)
    # 13 decline
    ls[3] = collect(2:12)
    # 15-18 decline
    ls[4] = collect(2:14)
    # 16 decline
    ls[5] = collect(2:15)
    # 20 decline
    ls[6] = collect(2:19)
    # no decline
    ls[7] = collect(4:20)
    # Then run multiple simulations
    a, J = prot_simulate_mult(ps,ai,Ni,Tmax)
    pyplot(dpi=200)
    m = L"^{-1}"
    pR = L"\phi_R"
    plot(xlabel="Energy acquisition rate ATP cell$(m) s$m",ylabel="ATP per cell")
    for i = 1:7
        # Tom used a log plot but I think this will obsurce too much at the moment
        scatter!(J[i,ls[i]],a[i,ls[i]],label="$(pR) = $(round(((i+1)/10)*0.55,digits=3))")
    end
    savefig("Output/ATPvsRate.png")
    plot(xlabel="Energy acquisition rate ATP cell$(m) s$m",ylabel="ATP per cell")
    for i = 1:7
        # Tom used a log plot but I think this will obsurce too much at the moment
        scatter!(log.(J[i,ls[i]]),log.(a[i,ls[i]]),label="$(pR) = $(round(((i+1)/10)*0.55,digits=3))")
    end
    savefig("Output/LogATPvsRate.png")
    return(nothing)
end

# function to make a plot of growth rate with changing ribosome fraction
function singpop_curv()
    println("Successfully compiled.")
    # Initialise parameter set
    ps = initialise_prot(false)
    # Pick an initial condition to optimise for (assume fixed environment)
    S = 100.0
    P = 10.0
    # Housekeeping fraction is fixed throughout
    ϕH = 0.45
    # Set number of fractions to consider
    N = 100
    # Preallocate output
    λ = zeros(N)
    ϕr = zeros(N)
    ϕ = zeros(3)
    # Fix housekeeping fraction
    ϕ[3] = ϕH
    for i = 1:N
        ϕ[1] = (1-ϕ[3])*(i-1)/(N-1)
        ϕ[2] = 1 - ϕ[1] - ϕ[3]
        # function to find maximum sustainable value for λ
        λ[i], _ = λ_max(S,P,ϕ,ps)
        ϕr[i] = ϕ[1]
    end
    pyplot(dpi=200)
    pR = L"\phi_R"
    ls = L"\lambda\;(s^{-1})"
    plot(ϕr,λ,label="",xlabel="Ribosome protein fraction, $(pR)",ylabel="Growth rate, $(ls)")
    savefig("Output/GrowthvsPR.png")
    # Make vector of ribosome fractions
    ϕR = collect(0:0.01:(1-ϕH))
    λ1 = zeros(length(ϕR))
    # Find growth rate with infinite energy for a fixed ribosome fraction
    a = 10e20 # energy taken to be approximatly infinite
    for i = 1:length(ϕR)
        λ1[i] = λs(a,ϕR[i],ps)
    end
    # Now need to implement second growth law, all energy goes directly to growth
    ϕP = collect((1-ϕH):-0.01:0)
    λ2 = zeros(length(ϕP))
    # Find growth rate assuming energy intake exactly balances energy use by translation
    for i = 1:length(ϕP)
        # Find amount of enzyme
        E = Eα(ϕP[i],ps)
        # Then use to find reaction rate
        q = qs(S::Float64,P::Float64,E,ps)
        J = ps.η*q
        λ2[i] = J/(ps.MC*ps.ρ)
    end
    plot(ϕR,λ1,label="Ribosome limited",xlabel="Ribosome protein fraction, $(pR)",ylabel="Growth rate, $(ls)")
    plot!(ϕR,λ2,label="Energy limited")
    savefig("Output/GrowthLaws.png")
    plot!(ϕr,λ,label="Actual")
    savefig("Output/CombGrowth.png")
    return(nothing)
end

# A testing function
function singpop_test()
    println("Successfully compiled.")
    # Initialise parameter set
    ps = initialise_prot(false)
    # Set amount of energy
    as = collect(1e7:1e7:1e9)
    ϕs = zeros(length(as))
    # Now find ev
    for i = 1:length(as)
        ϕs[i] = ϕ_R(as[i],ps)
    end
    plot(as,ϕs)
    savefig("Output/test.png")
    # Now find ribosome fraction vs a for a number of conditions
    return(nothing)
end

# Function to find range of responses to changing a
function singpop_range()
    println("Successfully compiled.")
    # Initialise parameter set
    ps = initialise_prot(true)
    S = 1e-2
    P = collect(0.0:1e-3:9e-3)
    ϕ = collect(0.01:0.01:0.55)
    θs = zeros(length(P))
    # Make vector of theta values
    for i = 1:length(P)
        θs[i] = θ(S,P[i],ps.T,ps.η,ps.r.ΔG0)
    end
    # Preallocate vectors for growth rate and energy concentration
    λ = zeros(length(P),length(ϕ))
    am = zeros(length(P),length(ϕ))
    # Now find λ and a values for each value of ϕ and P
    for j = 1:length(ϕ)
        for i = 1:length(P)
            λ[i,j], am[i,j] = λ_max(S,P[i],ϕ[j],ps)
        end
    end
    # The first ϕ value should have the highest range of a values
    ϕR = zeros(length(am[1,:]))
    # Find corresponding ϕR value for each step
    for i = 1:length(ϕR)
        ϕR[i] = ϕ_R(am[1,i],ps)
    end
    pyplot(dpi=200)
    plot(xlabel=L"\phi_R",ylabel="Energy")
    for i = 1:length(P)
        plot!(ϕ,am[i,:],label="θ = $(θs[i])")
    end
    plot!(ϕR,am[1,:],label="actual")
    savefig("Output/EnergyvsFraction.png")
    plot(xlabel=L"\phi_R",ylabel="Growth rate")
    for i = 1:length(P)
        plot!(ϕ,λ[i,:],label="θ = $(θs[i])")
    end
    savefig("Output/GrowthvsFraction.png")
end

# Simulation of a single population growing on a fixed initial amount of substrate
function singpop_fix()
    println("Successfully compiled.")
    # Simple test data set
    ai = 5.0 # initial energy level
    Ni = 10.0 # initial population
    Ci = 1.0
    # This is currently a paramter which I am fiddling
    Kγ = 1e8
    # Give omega a fairly arbitary value for now, would be expected to be of similar order to Kγ
    Ωs = [Kγ/100; Kγ/10; Kγ; 10*Kγ; 100*Kγ]
    # kc is another parameter that I should fiddle
    kc = 1.0
    for j = 1:length(Ωs)
        # Initialise parameter set
        ps = initialise_prot_fix(Kγ,Ωs[j],kc)
        # Choose simulation time
        Tmax = 250000.0
        # Then run simulation
        C, T = prot_simulate(ps,Tmax,ai,Ni,Ci)
        # Now calculate growth rates and proteome fractions
        λa = zeros(length(T))
        ϕR = zeros(length(T))
        for i = 1:length(T)
            ϕR[i] = ϕ_R(C[i,2],ps)
            λa[i] = λs(C[i,2],ϕR[i],ps)
        end
        # Do plotting
        pyplot(dpi=200)
        plot(T,C[:,1],xlabel="Time",label="",ylabel="Population")
        savefig("Output/FixPopvsTime$j.png")
        plot(T,C[:,2],xlabel="Time",label="",ylabel="Cell energy conc")
        savefig("Output/FixEnergyvsTime$j.png")
        plot(T,C[:,3:4],xlabel="Time",label=["Substrate" "Waste"],ylabel="Concentration")
        savefig("Output/FixMetabolitevsTime$j.png")
        s1 = L"s^{-1}"
        plot(T,λa,xlabel="Time",label="",ylabel="Growth rate $(s1)")
        savefig("Output/FixGrowthvsTime$j.png")
        plot(T,ϕR,xlabel="Time",label="",ylabel=L"\phi_R")
        plot!(T,C[:,5],label="")
        savefig("Output/FixFractionvsTime$j.png")
    end
    return(nothing)
end

@time singpop_fix()
