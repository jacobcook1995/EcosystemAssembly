# A script to run the dynamics for the full model.
export full_simulate, θ, θ_smooth, qs, test_full_simulate

export test_dynamics!, full_simulate_syn

# AT SOME POINT SOME OF THESE FUNCTIONS WILL NEED TO BE REWRITTEN

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

# version of θ function, that smmothes output by setting values > 1 to 1
function θ_smooth(S::Float64,P::Float64,T::Float64,η::Float64,ΔG0::Float64)
    # Catch perverse cases that sometimes arise
    if S <= 0.0
        θs = 1.0
    elseif P <= 0.0
        θs = 0.0
    else
        θs = Q(S,P)/Keq(T,η,ΔG0)
    end
    # In this case don't want to return θ values greater than 1
    return(min(θs,1.0))
end

# Extra function to find rate from thermodynamic inhibition, just useful for plotting
function qs(ps::Microbe,S::Float64,P::Float64,E::Float64,θs::Float64)
    q = ps.kc[1]*E*S*(1-θs)/(ps.KS[1] + S*(1+ps.kr[1]*θs))
    return(max(q,0.0))
end

# function to find the rate of substrate consumption by a particular reaction
function qs(S::Float64,P::Float64,E::Float64,i::Int64,ps::Microbe,T::Float64,r::Reaction)
    # To speed things I don't have a check here to ensure that r.ID matches ps.Reac[i]
    # This is something to check if I start getting errors
    θs = θ(S,P,T,ps.η[i],r.ΔG0)
    q = ps.kc[i]*E*S*(1-θs)/(ps.KS[i] + S*(1+ps.kr[i]*θs))
    # Ensure that negative value cannot be returned
    return(max(q,0.0))
end

# function to calculate the amount of a partcular enzyme a strain has
function Eα(ϕR::Float64,ps::Microbe,i::Int64)
    E = ps.MC*(1-ϕR-ps.ϕH)*ps.ϕP[i]/(ps.n[2])
    return(E)
end

# function to find the elongation rate γ
function γs(a::Float64,ps::Microbe)
    γ = ps.γm*a/(ps.Kγ+a)
    return(γ)
end

# function to find the growth rate λ
function λs(a::Float64,ϕR::Float64,ps::Microbe)
    # Find elongation rate
    γ = γs(a,ps)
    λ = (γ*ϕR*ps.Pb)/ps.n[1]
    return(λ)
end

# function to find ϕR based on the energy concentration
function ϕ_R(a::Float64,ps::Microbe)
    ϕ = (1-ps.ϕH)*a/(ps.KΩ + a)
    return(ϕ)
end

# function to implement the consumer resource dynamics
function full_dynamics!(dx::Array{Float64,1},x::Array{Float64,1},ms::Array{Microbe,1},ps::TOParameters,
                        rate::Array{Float64,2},t::Float64)
    # loop over the reactions to find reaction rate for each reaction for each strain
    for j = 1:ps.O
        # Find substrate and product for this reaction
        for i = 1:length(ms)
            # Check if microbe i performs reaction j
            if j ∈ ms[i].Reacs
                # Find index of this reaction in microbe
                k = findfirst(x->x==j,ms[i].Reacs)
                # Find amount of enzyme E
                E = Eα(x[2*length(ms)+ps.M+i],ms[i],k)
                # Then finally calculate reaction rate
                rate[i,j] = qs(x[length(ms)+ps.reacs[j].Rct],x[length(ms)+ps.reacs[j].Prd],E,k,ms[i],ps.T,ps.reacs[ms[i].Reacs[k]])
            else
                rate[i,j] = 0.0
            end
        end
    end
    # Now want to use the rate matrix in the consumer dynamics
    for i = 1:length(ms)
        # Check if strain is effectively extinct
        if x[i] <= 1e-10
            # If so x should be set to zero and should not change from that
            dx[i] = 0.0
            x[i] = 0.0
            # In this case the energy concentration should also be fixed to zero
            dx[length(ms)+ps.M+i] = 0.0
            x[length(ms)+ps.M+i] = 0.0
            # Corresponding proteome fraction also shouldn't shift
            dx[2*length(ms)+ps.M+i] = 0.0
        else
            # find growth rate for strains that aren't extinct
            λ = λs(x[length(ms)+ps.M+i],x[2*length(ms)+ps.M+i],ms[i])
            # (growth rate - death rate)*population
            dx[i] = (λ - ms[i].d)*x[i]
            # Now find optimal ribosome fraction
            ϕR = ϕ_R(x[length(ms)+ps.M+i],ms[i])
            # This introduces a time delay
            τ = ms[i].fd/λ
            # Then update actual ribosome fraction
            dx[2*length(ms)+ps.M+i] = (ϕR - x[2*length(ms)+ps.M+i])/τ
            # Energy intake is zero
            J = 0
            # Loop over all reactions to find energy gained by them
            for j = 1:ms[i].R
                J += ms[i].η[j]*rate[i,ms[i].Reacs[j]]
            end
            # Add energy intake and substract translation and dilution from the energy concentration
            dx[length(ms)+ps.M+i] = J - (ms[i].MC*ms[i].ρ + x[length(ms)+ps.M+i])*λ
        end
    end
    # Do basic resource dynamics
    for i = length(ms)+1:length(ms)+ps.M
        # fist add external supply of resource and decay
        dx[i] = ps.κ[i-length(ms)] - ps.δ[i-length(ms)]*x[i]
    end
    # Then loop over microbes
    for i = 1:length(ms)
        # Loop over reactions for specific microbe
        for j = 1:ms[i].R
            # Increase the product
            dx[length(ms)+ps.reacs[ms[i].Reacs[j]].Prd] += rate[i,ms[i].Reacs[j]]*x[i]/NA
            # and decrease the reactant
            dx[length(ms)+ps.reacs[ms[i].Reacs[j]].Rct] -= rate[i,ms[i].Reacs[j]]*x[i]/NA
        end
    end
    # Final step to correct for any concentrations that have dropped below threshold (1e-15)
    for i = length(ms)+1:length(ms)+ps.M
        # If the rate of change is above a threshold (1e-20) they are not altered
        if x[i] < 1e-15 && dx[i] <= 1e-20
            x[i] = 0.0
            dx[i] = 0.0
        end
    end
    # Any ATP numbers that have gone below 0.33 should be removed
    for i = (length(ms)+ps.M+1):(2*length(ms)+ps.M)
        if x[i] < 0.33
            x[i] = 0.0
            dx[i] = 0.0
        end
    end
    return(dx)
end

# NEED TO ADAPT THIS TO USE CALLBACKS + NEED TO ADAPT THIS TO PROPERLY USE NEW ps FUNCTION
# Simulation code to run one instatnce of the simulation with a user defined starting condition
# ps is parameter set, Tmax is the time to integrate to, pop, conc, as and ϕs are the intial conditions
# mpl is a pool of microbes
function full_simulate(ps::TOParameters,Tmax::Float64,pop::Float64,conc::Float64,as::Float64,ϕs::Float64,
                        mpl::Array{Microbe,1})
    # THIS PROCEDURE IS NECESSARY FOR TESTING, NEEDS TO BE CHANGED LATER, PARTICUALRLY TO AVOID DUPLICATES
    # Preallocate inital vector of microbes
    ms = Array{Microbe,1}(undef,10)
    # Randomly choose them from the pool
    for i = 1:length(ms)
        r = rand(1:length(mpl))
        ms[i] = mpl[r]
    end
    # Preallocate memory
    rate = zeros(length(ms),ps.O)
    # Now substitute preallocated memory in
    dyns!(dx,x,ms,t) = full_dynamics!(dx,x,ms,ps,rate,t)
    # Find time span for this step
    tspan = (0,Tmax)
    # Make
    pops = pop*ones(length(ms))
    concs = conc*ones(ps.M)
    ass = as*ones(length(ms))
    ϕss = ϕs*ones(length(ms))
    x0 = [pops;concs;ass;ϕss]
    # Then setup and solve the problem
    prob = ODEProblem(dyns!,x0,tspan,ms)
    # Still generates problems, not sure if I have to change a solver option or what
    sol = DifferentialEquations.solve(prob)
    return(sol',sol.t)
end
