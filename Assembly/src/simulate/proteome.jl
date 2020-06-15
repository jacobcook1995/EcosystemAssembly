# A script to test the protein fraction model for a single species and a single
# substrate. This script is for the initial implementation and testing of the
# proteome fraction model. Once I'm satisfied with it I will incorperate it into
# the main model

export prot_simulate, λs, optimise_ϕ, prot_simulate_mult, λ_max, Eα, qs, γs, ϕ_R, ϕ_R_2

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
function prot_simulate(ps::ProtParameters,Tmax::Float64,ai::Float64,Ni::Float64,Si::Float64=0.0)
    # Initialise vectors of concentrations and populations
    pop = Ni*ones(1)
    apop = ai*ones(1)
    conc = [Si,0.0,0.1] # No product to begin with, substrate can be set
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

# function to find maximum sustainable value for λ
function λ_max(S::Float64,P::Float64,ϕ::Float64,ps::ProtParameters)
    # If ribosome fraction is zero no growth is possible
    if ϕ == 0.0
        return(0.0,0.0)
    end
    # Calculate amount of enzyme
    E = Eα(1-ϕ-ps.ϕH,ps)
    # Then use rate to find J
    rate = qs(S,P,E,ps)
    # All energy comes from this reaction
    J = ps.η*rate
    # Define a as a symbol
    a = symbols("a")
    # Make full expression for dadt
    dadt = J - ((ps.γm*a)/(ps.Kγ+a))*(ϕ*ps.Pb/ps.n[1])*(ps.ρ*ps.MC+a)
    # Now solve this to find maximising a value
    af = convert(Float64,nsolve(dadt,1.0))
    # Find corresponding maximum rate
    λ = λs(af,ϕ,ps)
    return(λ,af)
end

# Function to find optimal ϕ for a given initial condition
function optimise_ϕ(S::Float64,P::Float64,ps::ProtParameters)
    # Preallocate output
    ϕm = zeros(3)
    # Seperate housekeeping fraction
    ϕT = 1 - ps.ϕH
    ϕm[3] = ps.ϕH
    # Make initial vector of ϕ
    ϕm[1] = (ϕT)/2
    ϕm[2] = 1 - ϕm[1] - ϕm[3]
    # Now find λ for this initial vector
    λm, af = λ_max(S,P,ϕm,ps)
    # loop until final value is found
    fine = false
    # Make new left and right phi functions
    ϕl = zeros(3)
    ϕr = zeros(3)
    ϕl[3] = ps.ϕH
    ϕr[3] = ps.ϕH
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
        λl, al = λ_max(S,P,ϕl,ps)
        λr, ar = λ_max(S,P,ϕr,ps)
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
            af = al
        else # otherwise go right
            ϕm = copy(ϕr)
            λm = λr
            af = ar
        end
    end
    return(ϕm,λm,af)
end

# function to simulate same species and return money
function prot_simulate_mult(ps::ProtParameters,ai::Float64,Ni::Float64,Ci::Float64,Tmax::Float64)
    # Initialise vectors of concentrations and populations
    pop = Ni*ones(1)
    apop = ai*ones(1)
    conc = [Ci; 0.0] # Large initial deposit of substrate
    ϕi = 0.1*ones(1) # Start with low ribosome fraction
    # Now sub the parameters in
    p_dyns!(dx,x,pa,t) = p_dynamics!(dx,x,pa,ps,t)
    # Make simulation time span
    tspan = (0,Tmax)
    x0 = [pop;apop;conc;ϕi]
    # parameter set kept empty
    pa = Array{Int64,1}(undef,0)
    # Then setup and solve the problem
    println("Simulation $(1) started.")
    prob = ODEProblem(p_dyns!,x0,tspan,pa)
    sol = DifferentialEquations.solve(prob)
    println("Final population $(sol'[end,1])")
    # Find position of the peak
    pk, tp = findmax(sol'[:,2])
    # find position where substrate has decayed below a certain point
    te = findlast(x->x>=0.9,sol'[:,3])
    # Preallocate output
    a = zeros(te+1-tp)
    J = zeros(te+1-tp)
    # Output all points up to final time point
    for j = tp:te
        a[j+1-tp] = sol'[j,2]
        E = Eα(1-sol'[j,5]-ps.ϕH,ps)
        J[j+1-tp] = ps.η*qs(sol'[j,3],sol'[j,4],E,ps)
    end
    return(a,J)
end
