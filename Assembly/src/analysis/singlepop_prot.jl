# A script to analyse the proteome model for a single population.
using Assembly
using Plots
using LaTeXStrings
using LsqFit
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
    C, T = prot_simulate(ps, Tmax, ai, Ni)
    # Now calculate growth rates and proteome fractions
    λa = zeros(length(T))
    ϕR = zeros(length(T))
    for i in 1:length(T)
        ϕR[i] = ϕ_R(C[i, 2], ps)
        λa[i] = λs(C[i, 2], ϕR[i], ps)
    end
    # Do plotting
    pyplot(dpi = 200, bg = :transparent, fg = :black)
    plot(T, C[:, 1], xlabel = "Time", label = "", ylabel = "Population", yaxis = :log)
    savefig("Output/PopvsTime.png")
    plot(T, C[:, 2], xlabel = "Time", label = "", ylabel = "Cell energy conc")
    savefig("Output/EnergyvsTime.png")
    plot(T, C[:, 3:4], xlabel = "Time", label = ["Substrate" "Waste"],
         ylabel = "Concentration")
    savefig("Output/MetabolitevsTime.png")
    s1 = L"s^{-1}"
    plot(T, λa, xlabel = "Time", label = "", ylabel = "Growth rate $(s1)")
    savefig("Output/GrowthvsTime.png")
    plot(T, ϕR, xlabel = "Time", label = "", ylabel = L"\phi_R")
    plot!(T, C[:, 5], label = "")
    savefig("Output/FractionvsTime.png")
    return (nothing)
end

# Same as above function but for a batch culture
function singpop_batch()
    println("Successfully compiled.")
    # Simple test data set
    ai = 5e6 # initial energy level
    Ni = 100.0 # initial population
    Si = 1.0
    # This is currently a parameter which I am fiddling
    Kγ = 5e8
    # Give omega a fairly arbitrary value for now, would be expected to be of similar order to Kγ
    KΩ = 2.0e9
    # kc is another parameter that I should fiddle
    kc = 1.0
    # Initialise parameter set
    ps = initialise_prot_fix(Kγ, KΩ, kc)
    # Choose simulation time
    Tmax = 250000.0
    # Then run simulation
    ϕi = 0.1 # Low initial ϕR
    C, T = prot_simulate(ps, Tmax, ai, Ni, Si, ϕi)
    # Now calculate growth rates and proteome fractions
    λa = zeros(length(T))
    ϕR = zeros(length(T))
    for i in 1:length(T)
        ϕR[i] = ϕ_R(C[i, 2], ps)
        λa[i] = λs(C[i, 2], ϕR[i], ps)
    end
    # Do plotting
    pyplot(dpi = 200, bg = :transparent, fg = :black)
    plot(T, C[:, 1], xlabel = "Time", label = "", ylabel = "Population", yaxis = :log)
    savefig("Output/PopvsTime.png")
    plot(T, C[:, 2], xlabel = "Time", label = "", ylabel = "ATP per cell")
    savefig("Output/EnergyvsTime.png")
    plot(T, C[:, 3:4], xlabel = "Time", label = ["Substrate" "Waste"],
         ylabel = "Concentration")
    savefig("Output/MetabolitevsTime.png")
    s1 = L"s^{-1}"
    plot(T, λa, xlabel = "Time", label = "", ylabel = "Growth rate $(s1)")
    savefig("Output/GrowthvsTime.png")
    pR = L"\phi_R"
    plot(T, ϕR, xlabel = "Time", label = "Optimal $(pR)", ylabel = L"\phi_R")
    plot!(T, C[:, 5], label = "Actual $(pR)")
    savefig("Output/FractionvsTime.png")
    pRi = L"\phi^i_R"
    plot(T, C[:, 2], xlabel = "Time", label = "$(pRi) = $(ϕi)", ylabel = "ATP per cell")
    savefig("Output/NotedCase.png")
    return (nothing)
end

function growth_laws()
    println("Successfully compiled.")
    # Simple test data set
    ai = 5.0 # initial energy level
    Ni = 100.0 # initial population
    Si = 1.0
    Pi = 0.0
    ϕi = 0.1
    # Choose simulation time
    Tmax = 5000000.0
    # Set this as a middling value of ΔG
    ΔG = -3e6 # Relatively small Gibbs free energy change
    # Max Elongation rate also taken from Bremer (1996), convert from minutes to seconds
    γm = 1260.0 / 60.0 # Change this one
    # Make vector of γm values
    γs = [
        γm / 100,
        γm / 50,
        γm / 20,
        γm / 10,
        γm / 7.5,
        γm / 5,
        γm / 4,
        γm / 3,
        γm / 2,
        γm / 1.5,
        γm,
    ]
    # Make vector to store final growth rates and fractions
    λ1 = zeros(length(γs))
    ϕ1 = zeros(length(γs))
    a1 = zeros(length(γs))
    # Setup plotting options
    pyplot(dpi = 200)
    pR = L"\phi_R"
    p1 = plot(xaxis = "Growth rate, λ", yaxis = "Ribosome fraction, $(pR)")
    for i in eachindex(γs)
        ps = initialise_prot_gl(γs[i], ΔG)
        C, T = prot_simulate_fix(ps, Tmax, ai, Ni, Si, Pi, ϕi)
        # Now calculate growth rates and proteome fractions
        λa = zeros(length(T))
        ϕR = zeros(length(T))
        for j in 1:length(T)
            ϕR[j] = ϕ_R(C[j, 2], ps)
            λa[j] = λs(C[j, 2], ϕR[j], ps)
        end
        # Save final λ and ϕR values
        λ1[i] = λa[end]
        ϕ1[i] = ϕR[end]
        a1[i] = C[end, 2]
    end
    # Now make set of ΔG values
    ΔGs = [
        ΔG,
        ΔG / 1.5,
        ΔG / 2,
        ΔG / 3,
        ΔG / 4,
        ΔG / 5,
        ΔG / 7.5,
        ΔG / 10,
        ΔG / 20,
        ΔG / 50,
        ΔG / 100,
    ]
    # Make vector to store final growth rates and fractions
    λ2 = zeros(length(ΔGs))
    ϕ2 = zeros(length(ΔGs))
    a2 = zeros(length(ΔGs))
    for i in eachindex(ΔGs)
        ps = initialise_prot_gl(γm, ΔGs[i])
        C, T = prot_simulate_fix(ps, Tmax, ai, Ni, Si, Pi, ϕi)
        # Now calculate growth rates and proteome fractions
        λa = zeros(length(T))
        ϕR = zeros(length(T))
        for j in 1:length(T)
            ϕR[j] = ϕ_R(C[j, 2], ps)
            λa[j] = λs(C[j, 2], ϕR[j], ps)
        end
        # Save final λ and ϕR values
        λ2[i] = λa[end]
        ϕ2[i] = ϕR[end]
        a2[i] = C[end, 2]
    end
    # Now want to do a least squares fit for both sets of data
    @. model(x, p) = p[1] * x + p[2]
    p0 = [0.5, 0.5]
    fit1 = curve_fit(model, λ1[4:end], ϕ1[4:end], p0)
    pr1 = coef(fit1)
    fit2 = curve_fit(model, λ2[1:(end - 3)], ϕ2[1:(end - 3)], p0)
    pr2 = coef(fit2)
    # plot both lines on the graph
    λ1s = [0.0; λ1]
    λ2s = [0.0; λ2]
    plot!(p1, λ1s, pr1[1] * λ1s .+ pr1[2], label = "", color = 1)
    plot!(p1, λ2s, pr2[1] * λ2s .+ pr2[2], label = "", color = 2)
    # Add arrows indicating direction of change
    l = 1e-4 # Way too large
    quiver!(p1, [λ1[end - 2]], [ϕ1[end - 2] + 0.03], quiver = ([-l], [-pr1[1] * l]),
            color = 3)
    quiver!(p1, [λ2[6]], [ϕ2[6] - 0.03], quiver = ([l], [pr2[1] * l]), color = 4)
    # Position is where the annotation centres are
    pos1x = λ1[end - 2] - l / 2
    pos1y = ϕ1[end - 2] + 0.04 - pr1[1] * l / 2
    pos2x = λ2[6] + l / 2
    pos2y = ϕ2[6] - 0.04 + pr2[1] * l / 2
    # Calculate rotations in degrees
    r1 = -12.5 # Needs to be -ve
    r2 = 20 # Needs to be +ve
    # Then add the annotations
    annotate!(p1, pos1x, pos1y, text("Translational inhibition", 5, :green, rotation = r1))
    annotate!(p1, pos2x, pos2y, text("Nutrient quality", 5, :purple, rotation = r2))
    # Plot final values
    scatter!(p1, λ1, ϕ1, label = "")
    scatter!(p1, λ2, ϕ2, label = "")
    # Make parameter set so that it can be used for plotting expressions
    ps = initialise_prot_gl(γm, ΔG)
    # Now plot lines from theory
    ϕ = collect(0.0:0.05:0.35)
    println("Theoretical slope = $(ps.γm*ps.Pb/ps.n[1])")
    println("Actual slope = $(1/pr2[1])")
    plot!(p1, (ps.γm * ps.Pb / ps.n[1]) * ϕ, ϕ, label = "Pred")
    # Finally save the graph
    savefig(p1, "Output/GrowthLaws.png")
    # Now plot the energy concentration data
    p2 = scatter(λ2, a2, label = "Nut Quality")
    # scatter!(p2,λ1,a1,label="Trans inhib")
    savefig(p2, "Output/EnergyData.png")
    return (nothing)
end

# function to compare the two different η competition cases
function ηcomparison()
    println("Successfully compiled.")
    # Changing now to investigate η competition
    N = 10
    Tmax = 100.0
    mq = 1.0
    mK = 0.1
    mk = 10.0
    # Run alternative simulations of the two cases
    ps1 = initialise_η(N, mq, mK, mk)
    C1, T1 = inhib_simulate(ps1, Tmax)
    # One extra strain as consumer exists
    ps2 = initialise_η2(N + 1, mq, mK, mk)
    C2, T2 = inhib_simulate(ps2, Tmax)
    # Run plotting
    pyplot(dpi = 200)
    plot(T1, C1[:, 1], label = "η = $(round(ps1.mics[1].η[1],digits=3))") # WRONG LABEL
    for i in 2:(ps1.N)
        plot!(T1, C1[:, i], label = "η = $(round(ps1.mics[i].η[1],digits=3))")
    end
    savefig("Output/Popsη1.png")
    plot(T1, C1[:, (ps1.N + 1):(ps1.N + 2)], label = "")
    savefig("Output/Concsη1.png")
    # Plot second case as well
    plot(T2, C2[:, 1], label = "consumer", color = ps2.N)
    for i in 2:(ps2.N)
        plot!(T2, C2[:, i], label = "η = $(round(ps2.mics[i].η[1],digits=3))",
              color = i - 1)
    end
    savefig("Output/Popsη2.png")
    plot(T2, C2[:, (ps2.N + 1):(ps2.N + 3)], label = "")
    savefig("Output/Concsη2.png")
    return (nothing)
end

# function to test the effect of varying cell mass MC
function singpop_MC()
    println("Successfully compiled.")
    # Simple test data set
    ai = 5.0 # initial energy level
    Ni = 100.0 # initial population
    MCs = [10^7, 10^8, 10^9]
    # Set up plotting
    pyplot(dpi = 200)
    p1 = plot(xlabel = "Time", ylabel = "Population", yaxis = :log)
    p2 = plot(xlabel = "Time", ylabel = "ATP molecules per cell")
    p3 = plot(xlabel = "Time", ylabel = "Concentration")
    s1 = L"s^{-1}"
    p4 = plot(xlabel = "Time", ylabel = "Growth rate $(s1)")
    p5 = plot(xlabel = "Time", ylabel = L"\phi_R")
    p6 = plot(xlabel = "Time", ylabel = "ATP molecules per amino acid")
    p7 = plot(xlabel = "Time", ylabel = "Total biomass", yaxis = :log)
    for i in eachindex(MCs)
        # Initialise parameter set
        ps = initialise_prot_M(MCs[i])
        # Choose simulation time
        Tmax = 1000000.0
        # Then run simulation
        C, T = prot_simulate(ps, Tmax, ai, Ni * (10^8 / MCs[i]))
        # Now calculate growth rates and proteome fractions
        λa = zeros(length(T))
        ϕR = zeros(length(T))
        for i in 1:length(T)
            ϕR[i] = ϕ_R(C[i, 2], ps)
            λa[i] = λs(C[i, 2], ϕR[i], ps)
        end
        # Do plotting
        plot!(p1, T, C[:, 1], label = "M = $(MCs[i])")
        plot!(p2, T, C[:, 2], label = "M = $(MCs[i])")
        plot!(p3, T, C[:, 3:4], label = ["Substrate" "Waste"])
        s1 = L"s^{-1}"
        plot!(p4, T, λa, label = "M = $(MCs[i])")
        plot!(p5, T, ϕR, label = "")
        plot!(p5, T, C[:, 5], label = "")
        plot!(p6, T, C[:, 2] / ps.MC, label = "M = $(MCs[i])")
        plot!(p7, T, C[:, 1] * ps.MC, label = "M = $(MCs[i])")
    end
    # Save the figures
    savefig(p1, "Output/PopvsTime.png")
    savefig(p2, "Output/EnergyvsTime.png")
    savefig(p3, "Output/MetabolitevsTime.png")
    savefig(p4, "Output/GrowthvsTime.png")
    savefig(p5, "Output/FractionvsTime.png")
    savefig(p6, "Output/RescaleEnergy.png")
    savefig(p7, "Output/TotalBiomass.png")
end

# function to test the effect of varying temperature T
function singpop_T()
    println("Successfully compiled.")
    # Simple test data set
    ai = 5.0 # initial energy level
    Ni = 100.0 # initial population
    Ts = [273.15, 283.15, 293.15, 303.15, 313.15]
    # Set up plotting
    pyplot(dpi = 200)
    p1 = plot(xlabel = "Time", ylabel = "Population", yaxis = :log)
    p2 = plot(xlabel = "Time", ylabel = "Cell energy conc")
    p3 = plot(xlabel = "Time", ylabel = "Concentration")
    s1 = L"s^{-1}"
    p4 = plot(xlabel = "Time", ylabel = "Growth rate $(s1)")
    p5 = plot(xlabel = "Time", ylabel = L"\phi_R")
    for i in eachindex(Ts)
        # Initialise product inhibited parameter set
        ps = initialise_prot_T(Ts[i], true)
        # Choose simulation time
        Tmax = 1000000.0
        # Then run simulation
        C, T = prot_simulate(ps, Tmax, ai, Ni)
        # Now calculate growth rates and proteome fractions
        λa = zeros(length(T))
        ϕR = zeros(length(T))
        for i in 1:length(T)
            ϕR[i] = ϕ_R(C[i, 2], ps)
            λa[i] = λs(C[i, 2], ϕR[i], ps)
        end
        # Do plotting
        plot!(p1, T, C[:, 1], label = "T = $(Ts[i]) K")
        plot!(p2, T, C[:, 2], label = "T = $(Ts[i]) K")
        plot!(p3, T, C[:, 3:4], label = ["Substrate" "Waste"])
        s1 = L"s^{-1}"
        plot!(p4, T, λa, label = "T = $(Ts[i]) K")
        plot!(p5, T, ϕR, label = "")
        plot!(p5, T, C[:, 5], label = "T = $(Ts[i]) K")
    end
    # Save the figures
    savefig(p1, "Output/PopvsTime.png")
    savefig(p2, "Output/EnergyvsTime.png")
    savefig(p3, "Output/MetabolitevsTime.png")
    savefig(p4, "Output/GrowthvsTime.png")
    savefig(p5, "Output/FractionvsTime.png")
    # Redo the plotting for noihibited case
    p1 = plot(xlabel = "Time", ylabel = "Population", yaxis = :log)
    p2 = plot(xlabel = "Time", ylabel = "Cell energy conc")
    p3 = plot(xlabel = "Time", ylabel = "Concentration")
    s1 = L"s^{-1}"
    p4 = plot(xlabel = "Time", ylabel = "Growth rate $(s1)")
    p5 = plot(xlabel = "Time", ylabel = L"\phi_R")
    for i in eachindex(Ts)
        # Initialise product inhibited parameter set
        ps = initialise_prot_T(Ts[i], false)
        # Choose simulation time
        Tmax = 1000000.0
        # Then run simulation
        C, T = prot_simulate(ps, Tmax, ai, Ni)
        # Now calculate growth rates and proteome fractions
        λa = zeros(length(T))
        ϕR = zeros(length(T))
        for i in 1:length(T)
            ϕR[i] = ϕ_R(C[i, 2], ps)
            λa[i] = λs(C[i, 2], ϕR[i], ps)
        end
        # Do plotting
        plot!(p1, T, C[:, 1], label = "T = $(Ts[i]) K")
        plot!(p2, T, C[:, 2], label = "T = $(Ts[i]) K")
        plot!(p3, T, C[:, 3:4], label = ["Substrate" "Waste"])
        s1 = L"s^{-1}"
        plot!(p4, T, λa, label = "T = $(Ts[i]) K")
        plot!(p5, T, ϕR, label = "")
        plot!(p5, T, C[:, 5], label = "T = $(Ts[i]) K")
    end
    # Save the figures
    savefig(p1, "Output/NoInhibPopvsTime.png")
    savefig(p2, "Output/NoInhibEnergyvsTime.png")
    savefig(p3, "Output/NoInhibMetabolitevsTime.png")
    savefig(p4, "Output/NoInhibGrowthvsTime.png")
    savefig(p5, "Output/NoInhibFractionvsTime.png")
end

# function to investigate synthesis rate vs efficiency tradeoff
function ρ_tradeoff()
    println("Compiled.")
    # Simple test data set
    ai = 5.0 # initial energy level
    Ni = 100.0 # initial population
    # Substrate and product for the optimisation
    S = 5.5e-3
    P = 5.5e-4
    fac = [1 / 4, 1 / 2, 3 / 4, 1] # factors to reduce γmax and ρ by
    # Set up plotting
    pyplot(dpi = 200)
    p1 = plot(xlabel = "Time", ylabel = "Population", yaxis = :log)
    p2 = plot(xlabel = "Time", ylabel = "Cell energy conc")
    p3 = plot(xlabel = "Time", ylabel = "Concentration")
    s1 = L"s^{-1}"
    p4 = plot(xlabel = "Time", ylabel = "Growth rate $(s1)")
    p5 = plot(xlabel = "Time", ylabel = L"\phi_R")
    # Now loop over the data
    for i in eachindex(fac)
        # Initialise parameter set
        ps = initialise_ρ(fac[i])
        # Update parameter set to have optimal K_Ω
        ps = opt_KΩ(ps, S, P)
        # Choose simulation time
        Tmax = 1000000.0
        # Then run simulation
        C, T = prot_simulate(ps, Tmax, ai, Ni)
        # Now calculate growth rates and proteome fractions
        λa = zeros(length(T))
        ϕR = zeros(length(T))
        for i in 1:length(T)
            ϕR[i] = ϕ_R(C[i, 2], ps)
            λa[i] = λs(C[i, 2], ϕR[i], ps)
        end
        # Do plotting
        plot!(p1, T, C[:, 1], label = "factor = $(fac[i])")
        plot!(p2, T, C[:, 2], label = "factor = $(fac[i])")
        plot!(p3, T, C[:, 3:4], label = ["Substrate" "Waste"])
        s1 = L"s^{-1}"
        plot!(p4, T, λa, label = "factor = $(fac[i])")
        plot!(p5, T, ϕR, label = "")
        plot!(p5, T, C[:, 5], label = "factor = $(fac[i])")
    end
    # Save the figures
    savefig(p1, "Output/EffPopvsTime.png")
    savefig(p2, "Output/EffEnergyvsTime.png")
    savefig(p3, "Output/EffMetabolitevsTime.png")
    savefig(p4, "Output/EffGrowthvsTime.png")
    savefig(p5, "Output/EffFractionvsTime.png")
    return (nothing)
end

# function to find the optimal value for KΩ for a particular parameter set
function opt_KΩ(ps::ProtParameters, S::Float64, P::Float64)
    # Choose simulation parameters
    Tmax = 1000000.0
    ai = 5.0
    Ni = 100.0
    ϕi = 0.1
    # Choose initial value of KΩ
    KΩ = 5e8
    ps = initialise_prot_KΩ(ps, KΩ)
    # Find optimal value for ϕR for this condition
    ϕr, _, am = optimise_ϕ(S, P, ps)
    # Run the simulations
    C, T = prot_simulate_fix(ps, Tmax, ai, Ni, S, P, ϕi)
    # Find final ribosome fraction
    ϕR = C[end, 5]
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
                δΩ = δΩ / 2
            end
        end
        # Make the two parameter sets
        psu = initialise_prot_KΩ(ps, Ku)
        psd = initialise_prot_KΩ(ps, Kd)
        # For some reason printing stops an error with DifferentialEquations????
        # println(KΩ)
        # Calculate final ϕ values for up and down cases
        C, T = prot_simulate_fix(psu, Tmax, ai, Ni, S, P, ϕi)
        ϕu = C[end, 5]
        C, T = prot_simulate_fix(psd, Tmax, ai, Ni, S, P, ϕi)
        ϕd = C[end, 5]
        # First check if already at the optimum
        if abs(ϕR - ϕr[1]) <= abs(ϕu - ϕr[1]) && abs(ϕR - ϕr[1]) <= abs(ϕd - ϕr[1])
            δΩ = δΩ / 2
            # Stop if step size is small and optimum has been reached
            if δΩ < 1.0e5
                fnd = true
            end
        elseif abs(ϕd - ϕr[1]) > abs(ϕu - ϕr[1]) # go up
            ϕR = ϕu
            KΩ = Ku
        else # otherwise go down
            ϕR = ϕd
            KΩ = Kd
        end
    end
    # Make final parameter set
    ps = initialise_prot_KΩ(ps, KΩ)
    return (ps)
end

@time singpop()
