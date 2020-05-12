# A script to test the protein fraction model for a single species and a single
# substrate. This script is for the initial implementation and testing of the
# proteome fraction model. Once I'm satisfied with it I will incorperate it into
# the main model
using Assembly
using DifferentialEquations # Needed as this is a test script not included in Assembly
using Plots
import PyPlot

# function to find the thermodynamic term θ, for the case of 1 to 1 stochiometry
function θ(S::Float64,P::Float64,T::Float64,η::Float64,ΔG0::Float64)
    # Catch perverse cases that sometimes arise
    if S <= 0.0
        θs = 1.0
    elseif P <= 0.0
        θs = 0.0
    else
        θs = Q(S,P)/Keq(T,η,ΔG0)
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

# Function to find optimal ϕ for a given initial condition
function optimise_ϕ(S::Float64,P::Float64,ϕH::Float64)
    # Preallocate output
    ϕ = zeros(3)
    # Seperate housekeeping fraction
    ϕT = 1 - ϕH
    ϕ[3] = ϕH

    return(ϕ)
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
    η = 0.9*(-ΔG/ΔGATP)
    k = 1.0 # matches qm from previously
    ai = 5.0 # initial energy level
    ρ = 1e-7
    n = [7459,300,300]
    # Pick an initial condition to optimise for (assume fixed environment)
    S = 50.0
    P = 5.0
    # Housekeeping protein fraction is constant
    ϕH = 0.2
    ϕ = optimise_ϕ(S,P,ϕH)

    ϕ = [0.4,0.4,0.2] # Again this should shift
    E = 1.0*(ϕ[2]/0.4)


    return(nothing)
end

@time singpop_opt()
