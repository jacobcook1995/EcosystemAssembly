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
        λ[i], _ = λ_max(S,P,ϕ[1],ps)
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
    Si = 1.0
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
        C, T = prot_simulate(ps,Tmax,ai,Ni,Si)
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

# Run a single population multiple times and plot a scatter graph of the results
function singpop_scat()
    println("Successfully compiled.")
    # Simple test data set
    ai = 5.0 # initial energy level
    Ni = 10.0 # initial population
    Ci = 1.0
    # This is currently a paramter which I am fiddling
    Kγ = 1e8
    # Give omega a fairly arbitary value for now, would be expected to be of similar order to Kγ
    Ωs = [Kγ]
    # kc is another parameter that I should fiddle
    kc = 1.0
    # choose eta values
    ηs = [0.9,1.1,1.2,1.25]*7.2
    # Choose simulation time
    Tmax = 250000.0
    # Setup plotting
    pyplot(dpi=200)
    m = L"^{-1}"
    pR = L"\phi_R"
    # Name the two plots
    p1 = plot(xlabel="Energy acquisition rate ATP cell$(m) s$m",ylabel="ATP per cell")
    p2 = plot(xlabel="Energy acquisition rate ATP cell$(m) s$m",ylabel="ATP per cell")
    # Now loop over the omega values
    for i = 1:length(Ωs)
        for j = 1:length(ηs)
            # Now make the parameter set
            ps = initialise_prot_fix(Kγ,Ωs[i],kc,ηs[j])
            # Then run multiple simulations
            a, J = prot_simulate_mult(ps,ai,Ni,Ci,Tmax)
            # Tom used a log plot but I think this will obsurce too much at the moment
            scatter!(p1,J,a,label="η = $(ηs[j])")
            # Tom used a log plot but I think this will obsurce too much at the moment
            scatter!(p2,log10.(J),log10.(a),label="η = $(ηs[j])")
        end
    end
    # Save the two plots
    savefig(p1,"Output/ATPvsRate.png")
    savefig(p2,"Output/LogATPvsRate.png")
    return(nothing)
end

# function to plot comparisons of different KΩ values
function singpop_KΩ_plots()
    println("Successfully compiled.")
    # Initialise parameter set
    ps = initialise_prot(false)
    # Pick an initial condition to optimise for (assume fixed environment)
    S = 5.5e-3
    P = 5.5e-4
    # Find optimal value for ϕR for this condition
    ϕr, _, am = optimise_ϕ(S,P,ps)
    # Choose simulation parameters
    Tmax = 1000000.0
    ai = 5.0
    Ni = 100.0
    ϕi = 0.1
    # Run simulation for this case
    C, T = prot_simulate_fix(ps,Tmax,ai,Ni,S,P,ϕr[1],false)
    # Now do plotting
    pyplot(dpi=200)
    p1 = plot(T,log10.(C[:,1]),xlabel="Time",label="",ylabel="Population")
    p2 = plot(T,C[:,2],xlabel="Time",label="",ylabel="Cell energy conc")
    p3 = plot(T,C[:,5],label="")
    # Now test for a range of Ω values
    Ωs = [ps.Kγ/10,ps.Kγ/5,ps.Kγ,5*ps.Kγ,10*ps.Kγ]
    # Loop over Ω values
    for j = 1:length(Ωs)
        # Make parameter set
        ps = initialise_prot_KΩ(Ωs[j])
        # Run the simulations
        C, T = prot_simulate_fix(ps,Tmax,ai,Ni,S,P,ϕi)
        # Then do the plotting
        plot!(p1,T,log10.(C[:,1]),label="Ω = $(Ωs[j])")
        plot!(p2,T,C[:,2],label="Ω = $(Ωs[j])")
        plot!(p3,T,C[:,5],label="Ω = $(Ωs[j])")
    end
    # Save all of the plots
    savefig(p1,"Output/PopvsTime.png")
    savefig(p2,"Output/EnergyvsTime.png")
    savefig(p3,"Output/FractionvsTime.png")
    return(nothing)
end

# function to find the optimal value for KΩ
function singpop_KΩ_opt()
    println("Successfully compiled.")
    # Initialise parameter set
    ps = initialise_prot(false)
    # Pick an initial condition to optimise for (assume fixed environment)
    S = 5.5e-3*100
    P = 5.5e-4
    # Find optimal value for ϕR for this condition
    ϕr, _, am = optimise_ϕ(S,P,ps)
    # Choose simulation parameters
    Tmax = 1000000.0
    ai = 5.0
    Ni = 100.0
    ϕi = 0.1
    # Choose intial value of KΩ
    KΩ = 5e8
    ps = initialise_prot_KΩ(KΩ)
    # Run the simulations
    C, T = prot_simulate_fix(ps,Tmax,ai,Ni,S,P,ϕi)
    # Find final ribosome fraction
    ϕR = C[end,5]
    println("Original ϕ = $(ϕr[1])")
    # Choose initial step size
    δΩ = 1e8
    # Preallocate variables so that they are in scope
    Ku = 5e8
    Kd = 5e8
    # Set up loop to find optimal KΩ value
    fnd = false
    while fnd == false
        # Step to generate new ϕ values
        update = false
        while update == false
            # Update up and down values
            Ku = KΩ + δΩ
            Kd = KΩ - δΩ
            if Kd >= 0.0
                update = true
            else
                # Reduce size of δϕ if it has gone negative
                δΩ = δΩ/2
            end
        end
        # Make the two parameter sets
        psu = initialise_prot_KΩ(Ku)
        psd = initialise_prot_KΩ(Kd)
        # Calculate final ϕ values for up and down cases
        C, T = prot_simulate_fix(psu,Tmax,ai,Ni,S,P,ϕi)
        ϕu = C[end,5]
        C, T = prot_simulate_fix(psd,Tmax,ai,Ni,S,P,ϕi)
        ϕd = C[end,5]
        # First check if already at the optimum
        if abs(ϕR-ϕr[1]) <= abs(ϕu-ϕr[1]) && abs(ϕR-ϕr[1]) <= abs(ϕd-ϕr[1])
            δΩ = δΩ/2
            # Stop if step size is small and optimum has been reached
            if δΩ < 1.0e5
                fnd = true
            end
        elseif abs(ϕd-ϕr[1]) > abs(ϕu-ϕr[1]) # go up
            ϕR = ϕu
            KΩ = Ku
        else # otherwise go down
            ϕR = ϕd
            KΩ = Kd
        end
    end
    println("Finished:")
    println("KΩ = $(KΩ)")
    println("ϕR = $(ϕR)")
    return(nothing)
end

# function to run comparisons of simulations with varying KΩ values
function singpop_KΩ_comp()
    println("Successfully compiled.")
    # Simple test data set
    ai = 5.0 # initial energy level
    Ni = 10.0 # initial population
    Si = 1.0
    # This is currently a parameter which I am fiddling
    Kγ = 5e8
    # Give omega a fairly arbitary value for now, would be expected to be of similar order to Kγ
    KΩs = [5e8,1e9,1.5e9,2e9,2.5e9,3e9,4e9,5e9,10e9]
    # kc is another parameter that I should fiddle
    kc = 1.0
    # Setup plotting
    pyplot(dpi=200)
    p1 = plot(xlabel="Time",ylabel="Population")
    p2 = plot(xlabel="Time",ylabel="Cell energy conc")
    p3 = plot(xlabel="Time",ylabel="Concentration")
    p4 = plot(xlabel="Time",ylabel=L"\phi_R")
    # Choose time interval
    Tmax = 2500000.0
    # Now loop over all simulations
    for i = 1:length(KΩs)
        ps = initialise_prot_fix(Kγ,KΩs[i],kc)
        # Run the simulation
        C, T = prot_simulate(ps,Tmax,ai,Ni,Si)
        # Plot results
        plot!(p1,T,log10.(C[:,1]),label="KΩ = $(KΩs[i])")
        plot!(p2,T,C[:,2],label="KΩ = $(KΩs[i])")
        plot!(p3,T,C[:,3:4],label="")
        plot!(p4,T,C[:,5],label="KΩ = $(KΩs[i])")
    end
    # Now save the plots
    savefig(p1,"Output/BatchPopvsTime.png")
    savefig(p2,"Output/BatchEnergyvsTime.png")
    savefig(p3,"Output/BatchMetabolitevsTime.png")
    savefig(p4,"Output/BatchFractionvsTime.png")
    # Now do chemostat set
    p1 = plot(xlabel="Time",ylabel="Population")
    p2 = plot(xlabel="Time",ylabel="Cell energy conc")
    p3 = plot(xlabel="Time",ylabel="Concentration")
    p4 = plot(xlabel="Time",ylabel=L"\phi_R")
    p5 = plot(xlabel="Time",ylabel="λ")
    p6 = plot(xlabel="Time",ylabel="J")
    p7 = plot(xlabel="Simulation number",ylabel="R*")
    p8 = plot(xlabel="Simulation number",ylabel="Time to max pop")
    # Choose time interval
    Tmax = 5000000.0
    # Now loop over all simulations
    for i = 1:length(KΩs)
        ps = initialise_prot_KΩ(KΩs[i],true)
        # Run the simulation
        C, T = prot_simulate(ps,Tmax,ai,Ni)
        # Calculate growth and energy aquistion rates
        λa = zeros(length(T))
        J = zeros(length(T))
        for i = 1:length(T)
            λa[i] = λs(C[i,2],C[i,5],ps)
            # Find amount of enzyme
            E = Eα(1-C[i,5]-ps.ϕH,ps)
            J[i] = ps.η*qs(C[i,3],C[i,4],E,ps)
        end
        # Plot results
        plot!(p1,T,log10.(C[:,1]),label="KΩ = $(KΩs[i])")
        plot!(p2,T,C[:,2],label="KΩ = $(KΩs[i])")
        plot!(p3,T,C[:,3],label="KΩ = $(KΩs[i])")
        plot!(p4,T,C[:,5],label="KΩ = $(KΩs[i])")
        plot!(p5,T,λa,label="KΩ = $(KΩs[i])")
        plot!(p6,T,J.*C[:,1],label="KΩ = $(KΩs[i])")
        scatter!(p7,[i],[C[end,3]],label="KΩ = $(KΩs[i])")
        # find maximum population
        mp, _ = findmax(C[:,1])
        # find first time point that is greater than 99% than this maximum
        tp = findfirst(x->x >= 0.99*mp,C[:,1])
        scatter!(p8,[i],[T[tp]],label="KΩ = $(KΩs[i])")
    end
    # Now save the plots
    savefig(p1,"Output/ChemoPopvsTime.png")
    savefig(p2,"Output/ChemoEnergyvsTime.png")
    savefig(p3,"Output/ChemoSubstratevsTime.png")
    savefig(p4,"Output/ChemoFractionvsTime.png")
    savefig(p5,"Output/ChemoGrowthvsTime.png")
    savefig(p6,"Output/ChemoGainvsTime.png")
    savefig(p7,"Output/ChemoEndSub.png")
    savefig(p8,"Output/ChemoEndTime.png")
    return(nothing)
end

@time singpop_KΩ_comp()
