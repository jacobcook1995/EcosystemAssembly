# A script to test the protein fraction model for a single species and a single
# substrate. This script is for the initial implementation and testing of the
# proteome fraction model. Once I'm satisfied with it I will incorperate it into
# the main model

export prot_simulate, λs, optimise_ϕ, prot_simulate_mult, λ_max, Eα, qs, γs, ϕ_R, ϕ_R_2, prot_simulate_fix

# function to find the rate of substrate consumption by a particular reaction
function qs(S::Float64,P::Float64,E::Float64,ps::ProtParameters)
    θs = θ(S,P,ps.T,ps.η,ps.r.ΔG0)
    q = ps.kc*E*S*(1-θs)/(ps.KS + S*(1+ps.kr*θs))
    # Ensure that negative value cannot be returned
    return(max(q,0.0))
end

# function to find the elongation rate γ
function γs(a::Float64,ps::ProtParameters)
    γ = ps.γm*a/(ps.Kγ+a)
    return(γ)
end

# function to find the growth rate λ
function λs(a::Float64,ϕR::Float64,ps::ProtParameters)
    # Find elongation rate
    γ = γs(a,ps)
    λ = (γ*ϕR*ps.Pb)/ps.n[1]
    return(λ)
end

# function to calculate the amount of a partcular enzyme a strain has
function Eα(ϕmα::Float64,ps::ProtParameters)
    E = ps.MC*ϕmα/(ps.n[2])
    return(E)
end

# function to find ϕR based on the energy concentration
function ϕ_R(a::Float64,ps::ProtParameters)
    ϕ = (1-ps.ϕH)*a/(ps.KΩ + a)
    return(ϕ)
end

# function to run the dynamics in the shifting proteome case
function p_dynamics!(dx::Array{Float64,1},x::Array{Float64,1},pa::Array{Int64,1},ps::ProtParameters,t::Float64)
    # Based on energy level find growth rate λ
    λ = λs(x[2],x[5],ps)
    # Then find amount of enzyme
    E = Eα(1-x[5]-ps.ϕH,ps)
    # Only one reaction so only one reaction rate
    rate = qs(x[3],x[4],E,ps)
    # All energy comes from this reaction
    J = ps.η*rate
    # Now update the stored energy
    dx[2] = J - (ps.MC*ps.ρ + x[2])*λ
    # Find optimal proteome fraction
    ϕR = ϕ_R(x[2],ps)
    # This introduces a time delay
    τ = ps.fd/λ
    dx[5] = (ϕR-x[5])/τ

    # Check if microbe is effectively dead
    if x[1] <= 1e-10
        # If so x should be set to zero and should not change from that
        dx[1] = 0.0
        x[1] = 0.0
    else
        # (growth rate - death rate)*population
        dx[1] = (λ - ps.d)*x[1]
    end
    # Do basic resource dynamics
    for i = 3:4
        # fist add external supply of resource and decay
        dx[i] = ps.κ[i-2] - ps.δ[i-2]*x[i]
    end
    # Only one reaction so consumption dynamics are simplified
    # Increase the product
    dx[4] += rate*x[1]/NA
    # and decrease the reactant
    dx[3] -= rate*x[1]/NA
    # Final step to correct for any concentrations that have gone negative
    for i = 2:4
        if x[i] < 0.0
            x[i] = 0.0
            dx[i] = 0.0
        end
    end
    return(dx)
end

# Simulation code to run one instatnce of the simulation
# ps is parameter set, Tmax is the time to integrate to
function prot_simulate(ps::ProtParameters,Tmax::Float64,ai::Float64,Ni::Float64,Si::Float64=0.0,ϕi::Float64=0.1)
    # Initialise vectors of concentrations and populations
    pop = Ni*ones(1)
    apop = ai*ones(1)
    conc = [Si,0.0,ϕi] # No product to begin with, substrate can be set
    # Now sub the parameters in
    p_dyns!(dx,x,pa,t) = p_dynamics!(dx,x,pa,ps,t)
    # Make simulation time span
    tspan = (0,Tmax)
    x0 = [pop;apop;conc]
    # Additional set of parameters so that it actually evaluates, empty for the moment
    pa = Array{Int64,1}(undef,0)
    # Then setup and solve the problem
    println("Simulation started.")
    prob = ODEProblem(p_dyns!,x0,tspan,pa)
    sol = DifferentialEquations.solve(prob)
    return(sol',sol.t)
end

# function to run the dynamics in the shifting proteome case with fixed concentrations
function p_dynamics_fix!(dx::Array{Float64,1},x::Array{Float64,1},pa::Array{Int64,1},ps::ProtParameters,t::Float64)
    # Based on energy level find growth rate λ
    λ = λs(x[2],x[5],ps)
    # Then find amount of enzyme
    E = Eα(1-x[5]-ps.ϕH,ps)
    # Only one reaction so only one reaction rate
    rate = qs(x[3],x[4],E,ps)
    # All energy comes from this reaction
    J = ps.η*rate
    # Now update the stored energy
    dx[2] = J - (ps.MC*ps.ρ + x[2])*λ
    # Find optimal proteome fraction
    ϕR = ϕ_R(x[2],ps)
    # This introduces a time delay
    τ = ps.fd/λ
    dx[5] = (ϕR-x[5])/τ

    # Check if microbe is effectively dead
    if x[1] <= 1e-10
        # If so x should be set to zero and should not change from that
        dx[1] = 0.0
        x[1] = 0.0
    else
        # (growth rate - death rate)*population
        dx[1] = (λ - ps.d)*x[1]
    end
    return(dx)
end

# function to run the dynamics in the shifting proteome case with fixed concentrations and no adaption
function p_dynamics_fix_na!(dx::Array{Float64,1},x::Array{Float64,1},pa::Array{Int64,1},ps::ProtParameters,t::Float64)
    # Based on energy level find growth rate λ
    λ = λs(x[2],x[5],ps)
    # Then find amount of enzyme
    E = Eα(1-x[5]-ps.ϕH,ps)
    # Only one reaction so only one reaction rate
    rate = qs(x[3],x[4],E,ps)
    # All energy comes from this reaction
    J = ps.η*rate
    # Now update the stored energy
    dx[2] = J - (ps.MC*ps.ρ + x[2])*λ

    # Check if microbe is effectively dead
    if x[1] <= 1e-10
        # If so x should be set to zero and should not change from that
        dx[1] = 0.0
        x[1] = 0.0
    else
        # (growth rate - death rate)*population
        dx[1] = (λ - ps.d)*x[1]
    end
    return(dx)
end

# Simulation code to run one instatnce of the simulation
# ps is parameter set, Tmax is the time to integrate to
function prot_simulate_fix(ps::ProtParameters,Tmax::Float64,ai::Float64,Ni::Float64,
                            Si::Float64,Pi::Float64,ϕi::Float64,upa::Bool=true)
    # Initialise vectors of concentrations and populations
    pop = Ni*ones(1)
    apop = ai*ones(1)
    conc = [Si,Pi,ϕi] # No product to begin with, substrate can be set
    if upa == true
        # Now sub the parameters in
        p_dyns1!(dx,x,pa,t) = p_dynamics_fix!(dx,x,pa,ps,t)
    else
        # MAYBE THIS IS CAUSING WARNINGS????
        p_dyns2!(dx,x,pa,t) = p_dynamics_fix_na!(dx,x,pa,ps,t)
    end
    # Make simulation time span
    tspan = (0,Tmax)
    x0 = [pop;apop;conc]
    # Additional set of parameters so that it actually evaluates, empty for the moment
    pa = Array{Int64,1}(undef,0)
    # Then setup and solve the problem
    # This splitting is needed to fix compilation
    if upa == true
        prob = ODEProblem(p_dyns1!,x0,tspan,pa)
    else
        prob = ODEProblem(p_dyns2!,x0,tspan,pa)
    end
    sol = DifferentialEquations.solve(prob)
    return(sol',sol.t)
end
