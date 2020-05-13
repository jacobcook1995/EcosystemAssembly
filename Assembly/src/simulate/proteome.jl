# A script to test the protein fraction model for a single species and a single
# substrate. This script is for the initial implementation and testing of the
# proteome fraction model. Once I'm satisfied with it I will incorperate it into
# the main model
using Assembly
using DifferentialEquations # Needed as this is a test script not included in Assembly
using SymPy # Again needed as this is a test script
using Plots
using LaTeXStrings
import PyPlot

# function to find the thermodynamic term θ, for the case of 1 to 1 stochiometry
function θ(S::Float64,P::Float64,T::Float64,η::Float64,ΔG0::Float64)
    # Catch perverse cases that sometimes arise
    if S <= 0.0
        θs = 1.0
    elseif P <= 0.0
        θs = 0.0
    else
        θs = Assembly.Q(S,P)/Keq(T,η,ΔG0) # NOT A FAN OF THIS NOTATION SHOULD OVERWRITE AT SOMEPOINT
    end
    # θ can be greater than 1, this does not have any impact as q cannot be negative
    return(θs)
end

# function to find the rate of substrate consumption by a particular reaction
function qs(S::Float64,P::Float64,T::Float64,η::Float64,ΔG0::Float64,k::Float64,E::Float64,KS::Float64,kr::Float64)
    θs = θ(S,P,T,η,ΔG0)
    q = k*E*S*(1-θs)/(KS + S*(1+kr*θs))
    # Ensure that negative value cannot be returned
    return(max(q,0.0))
end

# function to find the protein synthesis rate for a given protein fraction
function νx(a::Float64,n::Array{Int64,1},ϕ::Array{Float64,1},x::Int64)
    # HARD CODING THIS FOR NOW
    M = 10e8
    Pb = 0.5
    # Find elongation rate
    γ = γs(a)
    ν = (γ/n[x])*(M*ϕ[1]/n[1])*ϕ[x]*Pb
    return(ν)
end

# function to find the elongation rate γ
function γs(a::Float64)
    # HARD CODING THIS FOR NOW
    γm = 1260
    Kγ = 7
    γ = γm*a/(Kγ+a)
    return(γ)
end

# function to find the growth rate λ
function λs(a::Float64,nR::Int64,ϕR::Float64)
    # HARD CODING THIS FOR NOW
    Pb = 0.5
    # Find elongation rate
    γ = γs(a)
    λ = (γ*ϕR*Pb)/nR
    return(λ)
end

# function to run the dynamics in the shifting proteome case
function p_dynamics!(dx::Array{Float64,1},x::Array{Float64,1},ps::Array{Float64,1},T::Float64,d::Float64,
                    κ::Array{Float64,1},δ::Array{Float64,1},KS::Float64,kr::Float64,r::Reaction,η::Float64,
                    k::Float64,E::Float64,ϕ::Array{Float64,1},ρ::Float64,n::Array{Int64,1},t::Float64)
    # Only one reaction so only one reaction rate
    rate = qs(x[3],x[4],T,η,r.ΔG0,k,E,KS,kr)
    # All energy comes from this reaction
    J = η*rate
    T = 0
    for i = 1:length(ϕ)
        # LEAVE FOR NOW BUT THIS STEP IS ACTUALLLY SUPERFLUOUS
        T += ρ*n[i]*νx(x[2],n,ϕ,i)
    end
    # Then need to use this rate to find λ
    λ = λs(x[2],n[1],ϕ[1])
    # Now update the stored energy
    dx[2] = J - T - λ*x[2]

    # Check if microbe is effectively dead
    if x[1] <= 1e-10
        # If so x should be set to zero and should not change from that
        dx[1] = 0.0
        x[1] = 0.0
    else
        # (growth rate - death rate)*population
        dx[1] = (λ - d)*x[1]
    end
    # Do basic resource dynamics
    for i = 3:4
        # fist add external supply of resource and decay
        dx[i] = κ[i-2] - δ[i-2]*x[i]
    end
    # Only one reaction so consumption dynamics are simplified
    # Increase the product
    dx[4] += rate*x[1]
    # and decrease the reactant
    dx[3] -= rate*x[1]
    # Final step to correct for any concentrations that have gone negative
    for i = 2:4
        if x[i] < 0.0
            x[i] = 0.0
            dx[i] = 0.0
        end
    end
    return(dx)
end

# Just want to set up and run a single population
function singpop()
    # Simple test data set
    T = 293.15 # Assume that temperature T is constant at 20°C
    d = 5e-3 # death rate
    κ = [100.0,0.0] # Metabolite supply rate
    δ = 1.0*ones(2) # Metabolite dilution rate
    KS = 0.1 # Saturation constant
    kr = 10.0 # Reversibility factor
    # Now make the reaction
    ΔG = -6e5 # Relatively small Gibbs free energy change
    r = make_Reaction(1,1,2,ΔG)
    η = 0.9*(-ΔG/ΔGATP)
    k = 1.0 # matches qm from previously
    E = 1.0 # FIXED FOR NOW BUT SHOULD BE ADJUSTED LATER
    ai = 5.0 # initial energy level
    ϕ = [0.4,0.4,0.2] # Again this should shift
    ρ = 1e-7
    n = [7459,300,300]
    # Initialise vectors of concentrations and populations
    pop = 100*ones(1)
    apop = ai*ones(1)
    conc = zeros(2) # No chemical to begin with
    # Now sub the parameters in
    p_dyns!(dx,x,ps,t) = p_dynamics!(dx,x,ps,T,d,κ,δ,KS,kr,r,η,k,E,ϕ,ρ,n,t)
    # Not actually using parameters but need it to run
    ps = zeros(0)
    # Find time span for this step
    Tmax = 2500.0
    tspan = (0,Tmax)
    x0 = [pop;apop;conc]
    # Then setup and solve the problem
    println("Simulation started.")
    prob = ODEProblem(p_dyns!,x0,tspan,ps)
    sol = DifferentialEquations.solve(prob)
    # Do plotting
    pyplot(dpi=200)
    plot(sol.t,sol'[:,1],xlabel="Time",label="",ylabel="Population")
    savefig("Output/testPop.png")
    plot(sol.t,sol'[:,2],xlabel="Time",label="",ylabel="Cell energy conc")
    savefig("Output/testEng.png")
    plot(sol.t,sol'[:,3:4],xlabel="Time",label=["Substrate" "Waste"],ylabel="Concentration")
    savefig("Output/testCon.png")
    return(nothing)
end

# function to find maximum sustainable value for λ
function λ_max(S::Float64,P::Float64,ϕ::Array{Float64,1},r::Reaction,T::Float64,η::Float64,
                k::Float64,KS::Float64,kr::Float64,nR::Int64,ρ::Float64)
    # HARD CODING THIS FOR NOW
    γm = 1260
    Kγ = 7
    Pb = 0.5
    M = 10e8
    # first use ϕ to find E
    E = 1.0*(ϕ[2]/0.4) # THIS NEEDS TO BE BETTER WRITTEN LATER, AS A FUNCTION
    # Then use rate to find J
    rate = qs(S,P,T,η,r.ΔG0,k,E,KS,kr)
    # All energy comes from this reaction
    J = η*rate
    # Define a as a symbol
    a = symbols("a")
    # Make full expression for dadt
    dadt = J - ((γm*a)/(Kγ+a))*(ϕ[1]*Pb/nR)*(ρ*M+a)
    # Now solve this to find maximising a value
    af = convert(Float64,nsolve(dadt,1.0))
    # Find corresponding maximum rate
    λ = λs(af,nR,ϕ[1])
    return(λ)
end

# Function to find optimal ϕ for a given initial condition
function optimise_ϕ(S::Float64,P::Float64,ϕH::Float64,r::Reaction,T::Float64,η::Float64,
                    k::Float64,KS::Float64,kr::Float64,nR::Int64,ρ::Float64)
    # Preallocate output
    ϕm = zeros(3)
    # Seperate housekeeping fraction
    ϕT = 1 - ϕH
    ϕm[3] = ϕH
    # Make initial vector of ϕ
    ϕm[1] = (ϕT)/2
    ϕm[2] = 1 - ϕm[1] - ϕm[3]
    # Now find λ for this initial vector
    λm = λ_max(S,P,ϕm,r,T,η,k,KS,kr,nR,ρ)
    # loop until final value is found
    fine = false
    # Make new left and right phi functions
    ϕl = zeros(3)
    ϕr = zeros(3)
    ϕl[3] = ϕH
    ϕr[3] = ϕH
    δϕ = 0.1
    while fine == false
        # Step to generate new ϕ values
        update = false
        while update == false
            # Update left and right functions
            ϕl[1] = ϕm[1] + δϕ
            ϕl[2] = 1 - ϕl[1] - ϕl[3]
            ϕr[2] = ϕm[2] + δϕ
            ϕr[1] = 1 - ϕr[2] - ϕr[3]
            if ϕl[2] >= 0.0 && ϕr[1] >= 0.0
                update = true
            else
                # Reduce size of δϕ if it has gone negative
                δϕ = δϕ/2
            end
        end
        # Calculate left and right growth rates
        λl = λ_max(S,P,ϕl,r,T,η,k,KS,kr,nR,ρ)
        λr = λ_max(S,P,ϕr,r,T,η,k,KS,kr,nR,ρ)
        # First check if already at the optimum
        if λm >= λl && λm >= λr
            δϕ = δϕ/2
            # Stop if step size is small and optimum has been reached
            if δϕ < 1.0e-5
                fine = true
            end
        elseif λl > λr # go left
            ϕm = copy(ϕl)
            λm = λl
        else # otherwise go right
            ϕm = copy(ϕr)
            λm = λr
        end
    end
    return(ϕm,λm)
end

# Similar function as above except that it tries to find optimal values for ϕ_M
function singpop_opt()
    # Simple test data set
    T = 293.15 # Assume that temperature T is constant at 20°C
    d = 5e-3 # death rate
    κ = [100.0,0.0] # Metabolite supply rate
    δ = 1.0*ones(2) # Metabolite dilution rate
    KS = 0.1 # Saturation constant
    kr = 10.0 # Reversibility factor
    # Now make the reaction
    ΔG = -6e5 # Relatively small Gibbs free energy change
    r = make_Reaction(1,1,2,ΔG)
    η = 1.0*(-ΔG/ΔGATP)
    k = 1.0 # matches qm from previously
    ai = 5.0 # initial energy level
    ρ = 1e-7
    n = [7459,300,300]
    # Pick an initial condition to optimise for (assume fixed environment)
    S = 100.0
    P = collect(5.0:5.0:95.0)
    θs = zeros(length(P))
    # Make vector of theta values
    for i = 1:length(P)
        θs[i] = θ(S,P[i],T,η,ΔG)
    end
    # Now find optimal ribosome fractions and growth rates for each one
    ϕR = zeros(length(P))
    λs = zeros(length(P))
    # Housekeeping protein fraction is constant
    ϕH = 0.2
    for i = 1:length(P)
        ϕ, λs[i] = optimise_ϕ(S,P[i],ϕH,r,T,η,k,KS,kr,n[1],ρ)
        ϕR[i] = ϕ[1]
    end
    pyplot(dpi=200)
    plot(θs,ϕR,label="",xlabel=L"θ",ylabel=L"ϕ_R")
    savefig("Output/RibosomeFrac.png")
    plot(θs,λs,label="",xlabel=L"θ",ylabel="Optimal growth rate")
    savefig("Output/GrowthRate.png")
    return(nothing)
end

@time singpop_opt()
