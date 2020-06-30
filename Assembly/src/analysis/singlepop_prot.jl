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
    plot(T,C[:,1],xlabel="Time",label="",ylabel="Population",yaxis=:log)
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
    γm = 1260.0/60.0 # Change this one
    # Make vector of γm values
    γs = [γm/100,γm/50,γm/20,γm/10,γm/7.5,γm/5,γm/4,γm/3,γm/2,γm/1.5,γm]
    # Make vector to store final growth rates and fractions
    λ1 = zeros(length(γs))
    ϕ1 = zeros(length(γs))
    # Setup plotting options
    pyplot(dpi=200)
    pR = L"\phi_R"
    p1 = plot(xaxis="Growth rate, λ",yaxis="Ribosome fraction, $(pR)")
    for i = 1:length(γs)
        ps = initialise_prot_gl(γs[i],ΔG)
        C, T = prot_simulate_fix(ps,Tmax,ai,Ni,Si,Pi,ϕi)
        # Now calculate growth rates and proteome fractions
        λa = zeros(length(T))
        ϕR = zeros(length(T))
        for j = 1:length(T)
            ϕR[j] = ϕ_R(C[j,2],ps)
            λa[j] = λs(C[j,2],ϕR[j],ps)
        end
        # Save final λ and ϕR values
        λ1[i] = λa[end]
        ϕ1[i] = ϕR[end]
    end
    # Now make set of ΔG values
    ΔGs = [ΔG,ΔG/1.5,ΔG/2,ΔG/3,ΔG/4,ΔG/5,ΔG/7.5,ΔG/10,ΔG/20,ΔG/50,ΔG/100]
    # Make vector to store final growth rates and fractions
    λ2 = zeros(length(ΔGs))
    ϕ2 = zeros(length(ΔGs))
    for i = 1:length(ΔGs)
        ps = initialise_prot_gl(γm,ΔGs[i])
        C, T = prot_simulate_fix(ps,Tmax,ai,Ni,Si,Pi,ϕi)
        # Now calculate growth rates and proteome fractions
        λa = zeros(length(T))
        ϕR = zeros(length(T))
        for j = 1:length(T)
            ϕR[j] = ϕ_R(C[j,2],ps)
            λa[j] = λs(C[j,2],ϕR[j],ps)
        end
        # Save final λ and ϕR values
        λ2[i] = λa[end]
        ϕ2[i] = ϕR[end]
    end
    # Now want to do a least squares fit for both sets of data
    @. model(x, p) = p[1]*x + p[2]
    p0 = [0.5,0.5]
    fit1 = curve_fit(model,λ1[4:end],ϕ1[4:end],p0)
    pr1 = coef(fit1)
    fit2 = curve_fit(model,λ2[1:end-3],ϕ2[1:end-3],p0)
    pr2 = coef(fit2)
    # plot both lines on the graph
    λ1s = [0.0;λ1]
    λ2s = [0.0;λ2]
    plot!(p1,λ1s,pr1[1]*λ1s.+pr1[2],label="",color=1)
    plot!(p1,λ2s,pr2[1]*λ2s.+pr2[2],label="",color=2)
    # Add arrows indicating direction of change
    l = 1e-4 # Way too large
    quiver!(p1,[λ1[end-2]],[ϕ1[end-2]+0.03],quiver=([-l],[-pr1[1]*l]),color=3)
    quiver!(p1,[λ2[6]],[ϕ2[6]-0.03],quiver=([l],[pr2[1]*l]),color=4)
    # Position is where the annotation centres are
    pos1x = λ1[end-2] - l/2
    pos1y = ϕ1[end-2] + 0.04 - pr1[1]*l/2
    pos2x = λ2[6] + l/2
    pos2y = ϕ2[6] - 0.04 + pr2[1]*l/2
    # Calculate rotations in degrees
    r1 = -12.5 # Needs to be -ve
    r2 = 20 # Needs to be +ve
    # Then add the annotations
    annotate!(p1,pos1x,pos1y,text("Translational inhibition",5,:green,rotation=r1))
    annotate!(p1,pos2x,pos2y,text("Nutrient quality",5,:purple,rotation=r2))
    # Plot final values
    scatter!(p1,λ1,ϕ1,label="")
    scatter!(p1,λ2,ϕ2,label="")
    # Finally save the graph
    savefig(p1,"Output/GrowthLaws.png")
    return(nothing)
end

@time growth_laws()
