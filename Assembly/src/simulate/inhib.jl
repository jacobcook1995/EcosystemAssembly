# Script that attempts to implement an alternative version of the Marsland model that
# utilizes themodynamic end product inhibition
using Assembly

export inhib_simulate, test_inhib_simulate

# function to find the reaction quotient Q, in the case of 1 to 1 stochiometery
function Q(S::Float64,P::Float64)
    Q = P/S
    return(Q)
end

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
function qs(S::Float64,P::Float64,T::Float64,η::Float64,ΔG0::Float64,qm::Float64,KS::Float64,kr::Float64)
    θs = θ(S,P,T,η,ΔG0)
    q = qm*S*(1-θs)/(KS + S*(1+kr*θs))
    # Ensure that negative value cannot be returned
    return(max(q,0.0))
end


# function to implement the consumer resource dynamics
function dynamics!(dx::Array{Float64,1},x::Array{Float64,1},ps::InhibParameters,rate::Array{Float64,2},t::Float64)
    # loop over the reactions to find reaction rate for each reaction for each strain
    for j = 1:ps.O
        # Find substrate and product for this reaction
        for i = 1:ps.N
            # Check if microbe i performs reaction j
            if j ∈ ps.mics[i].Reacs
                # Find index of this reaction in microbe
                k = findfirst(x->x==j,ps.mics[i].Reacs)
                # Use k to select correct kinetic parameters
                rate[i,j] = qs(x[ps.N+ps.reacs[j].Rct],x[ps.N+ps.reacs[j].Prd],ps.T,ps.mics[i].η[k],ps.reacs[j].ΔG0,ps.mics[i].qm[k],ps.mics[i].KS[k],ps.mics[i].kr[k])
            else
                rate[i,j] = 0.0
            end
        end
    end
    # Now want to use the rate matrix in the consumer dynamics
    for i = 1:ps.N
        # Check if microbe is effectively dead
        if x[i] <= 1e-10
            # If so x should be set to zero and should not change from that
            dx[i] = 0.0
            x[i] = 0.0
        else
            # subtract maintenance
            dx[i] = -ps.mics[i].m
            # Add all the non zero contributions
            for j = 1:ps.mics[i].R
                dx[i] += ps.mics[i].η[j]*rate[i,ps.mics[i].Reacs[j]]
            end
            # multiply by population and proportionality constant
            dx[i] *= ps.mics[i].g*x[i]
        end
    end
    # Do basic resource dynamics
    for i = ps.N+1:ps.N+ps.M
        # fist add external supply of resource and decay
        dx[i] = ps.κ[i-ps.N] - ps.δ[i-ps.N]*x[i]
    end
    # Then loop over microbes
    for i = 1:ps.N
        # Loop over reactions for specific microbe
        for j = 1:ps.mics[i].R
            # Increase the product
            dx[ps.N+ps.reacs[ps.mics[i].Reacs[j]].Prd] += rate[i,ps.mics[i].Reacs[j]]*x[i]
            # and decrease the reactant
            dx[ps.N+ps.reacs[ps.mics[i].Reacs[j]].Rct] -= rate[i,ps.mics[i].Reacs[j]]*x[i]
        end
    end
    # Final step to correct for any concentrations that have gone negative
    for i = ps.N+1:ps.N+ps.M
        if x[i] < 0.0
            x[i] = 0.0
            dx[i] = 0.0
        end
    end
    return(dx)
end

# Simulation code to run one instatnce of the simulation
# ps is parameter set, Tmax is the time to integrate to
function inhib_simulate(ps::InhibParameters,Tmax::Float64)
    # Initialise vectors of concentrations and populations
    pop = ones(ps.N)
    conc = zeros(ps.M) # No chemical to begin with
    # Preallocate memory
    rate = zeros(ps.N,ps.O)
    # Now substitute preallocated memory in
    dyns!(dx,x,ps,t) = dynamics!(dx,x,ps,rate,t)
    # Find time span for this step
    tspan = (0,Tmax)
    x0 = [pop;conc]
    # Then setup and solve the problem
    println("Simulation started.")
    prob = ODEProblem(dyns!,x0,tspan,ps)
    sol = DifferentialEquations.solve(prob)
    return(sol',sol.t)
end

# Same dynamics function altered to do detailed testing
function test_dynamics!(dx::Array{Float64,1},x::Array{Float64,1},ps::InhibParameters,rate::Array{Float64,2},t::Float64)
    # loop over the reactions to find reaction rate for each reaction for each strain
    for j = 1:ps.O
        # Find substrate and product for this reaction
        for i = 1:ps.N
            # Check if microbe i performs reaction j
            if j ∈ ps.mics[i].Reacs
                # Find index of this reaction in microbe
                k = findfirst(x->x==j,ps.mics[i].Reacs)
                # Use k to select correct kinetic parameters
                rate[i,j] = qs(x[ps.N+ps.reacs[j].Rct],x[ps.N+ps.reacs[j].Prd],ps.T,ps.mics[i].η[k],ps.reacs[j].ΔG0,ps.mics[i].qm[k],ps.mics[i].KS[k],ps.mics[i].kr[k])
            else
                rate[i,j] = 0.0
            end
        end
    end
    # Now want to use the rate matrix in the consumer dynamics
    for i = 1:ps.N
        # Check if microbe is effectively dead
        if x[i] <= 1e-10
            # If so x should be set to zero and should not change from that
            dx[i] = 0.0
            x[i] = 0.0
        else
            # subtract maintenance
            dx[i] = -ps.mics[i].m
            # Add all the non zero contributions
            for j = 1:ps.mics[i].R
                dx[i] += ps.mics[i].η[j]*rate[i,ps.mics[i].Reacs[j]]
            end
            # multiply by population and proportionality constant
            dx[i] *= ps.mics[i].g*x[i]
        end
    end
    # Do basic resource dynamics
    for i = ps.N+1:ps.N+ps.M
        # fist add external supply of resource and decay
        dx[i] = ps.κ[i-ps.N] - ps.δ[i-ps.N]*x[i]
    end
    # Then loop over microbes
    for i = 1:ps.N
        # Loop over reactions for specific microbe
        for j = 1:ps.mics[i].R
            # Increase the product
            dx[ps.N+ps.reacs[ps.mics[i].Reacs[j]].Prd] += rate[i,ps.mics[i].Reacs[j]]*x[i]
            # and decrease the reactant
            dx[ps.N+ps.reacs[ps.mics[i].Reacs[j]].Rct] -= rate[i,ps.mics[i].Reacs[j]]*x[i]
        end
    end
    # Final step to correct for any concentrations that have gone negative
    for i = ps.N+1:ps.N+ps.M
        if x[i] < 0.0
            x[i] = 0.0
            dx[i] = 0.0
        end
    end
    return(dx)
end

# Same function as above but for detailed testing
function test_inhib_simulate(ps::InhibParameters,Tmax::Float64)
    # Initialise vectors of concentrations and populations
    pop = ones(ps.N)
    conc = zeros(ps.M) # No chemical to begin with
    # Preallocate memory
    rate = zeros(ps.N,ps.O)
    # Now substitute preallocated memory in
    dyns!(dx,x,ps,t) = test_dynamics!(dx,x,ps,rate,t)
    # Find time span for this step
    tspan = (0,Tmax)
    x0 = [pop;conc]
    # Then setup and solve the problem
    println("Test simulation started.")
    prob = ODEProblem(dyns!,x0,tspan,ps)
    sol = DifferentialEquations.solve(prob)
    println(sol.destats) # Also useful
    println(sol.retcode) # Useful
    return(sol',sol.t)
end
