# A script to run the dynamics for the full model.
export full_simulate, θ, θ_smooth, qs, test_full_simulate

export test_dynamics!, full_simulate_syn

# function to find the thermodynamic term θ, for the case of 1 to 1 stoichiometry
function θ(S::Float64, P::Float64, T::Float64, η::Float64, ΔG0::Float64)
    # Catch perverse cases that sometimes arise
    if S <= 0.0
        θs = 1.0
    elseif P <= 0.0
        θs = 0.0
    else
        θs = Q(S, P) / Keq(T, η, ΔG0)
    end
    # θ can be greater than 1, this does not have any impact as q cannot be negative
    return (θs)
end

# version of θ function, that smooths output by setting values > 1 to 1
function θ_smooth(S::Float64, P::Float64, T::Float64, η::Float64, ΔG0::Float64)
    # Catch perverse cases that sometimes arise
    if S <= 0.0
        θs = 1.0
    elseif P <= 0.0
        θs = 0.0
    else
        θs = Q(S, P) / Keq(T, η, ΔG0)
    end
    # In this case don't want to return θ values greater than 1
    return (min(θs, 1.0))
end

# Extra function to find rate from thermodynamic inhibition, just useful for plotting
function qs(ps::MicrobeP, S::Float64, P::Float64, E::Float64, θs::Float64)
    q = ps.kc[1] * E * S * (1 - θs) / (ps.KS[1] + S * (1 + ps.kr[1] * θs))
    return (max(q, 0.0))
end

# function to find the rate of substrate consumption by a particular reaction
function qs(S::Float64, P::Float64, E::Float64, i::Int64, ps::MicrobeP, T::Float64,
            r::Reaction)
    # To speed things I don't have a check here to ensure that r.ID matches ps.Reac[i]
    # This is something to check if I start getting errors
    θs = θ(S, P, T, ps.η[i], r.ΔG0)
    q = ps.kc[i] * E * S * (1 - θs) / (ps.KS[i] + S * (1 + ps.kr[i] * θs))
    # Ensure that negative value cannot be returned
    return (max(q, 0.0))
end

# function to calculate the amount of a particular enzyme a strain has
function Eα(ϕR::Float64, ps::MicrobeP, i::Int64)
    E = ps.MC * (1 - ϕR - ps.ϕH) * ps.ϕP[i] / (ps.n[2])
    return (E)
end

# function to find the elongation rate γ
function γs(a::Float64, ps::MicrobeP)
    γ = ps.γm * a / (ps.Kγ + a)
    return (γ)
end

# function to find the growth rate λ
function λs(a::Float64, ϕR::Float64, ps::MicrobeP)
    # Find elongation rate
    γ = γs(a, ps)
    λ = (γ * ϕR * ps.Pb) / ps.n[1]
    return (λ)
end

# function to find ϕR based on the energy concentration
function ϕ_R(a::Float64, ps::MicrobeP)
    ϕ = (1 - ps.ϕH) * a / (ps.KΩ + a)
    return (ϕ)
end

# function to implement the consumer resource dynamics
function full_dynamics!(dx::Array{Float64, 1}, x::Array{Float64, 1}, ps::FullParameters,
                        rate::Array{Float64, 2}, t::Float64)
    # loop over the reactions to find reaction rate for each reaction for each strain
    for j in 1:(ps.O)
        # Find substrate and product for this reaction
        for i in 1:(ps.N)
            # Check if microbe i performs reaction j
            if j ∈ ps.mics[i].Reacs
                # Find index of this reaction in microbe
                k = findfirst(x -> x == j, ps.mics[i].Reacs)
                # Find amount of enzyme E
                E = Eα(x[2 * ps.N + ps.M + i], ps.mics[i], k)
                # Then finally calculate reaction rate
                rate[i, j] = qs(x[ps.N + ps.reacs[j].Rct], x[ps.N + ps.reacs[j].Prd], E, k,
                                ps.mics[i], ps.T, ps.reacs[ps.mics[i].Reacs[k]])
            else
                rate[i, j] = 0.0
            end
        end
    end
    # Now want to use the rate matrix in the consumer dynamics
    for i in 1:(ps.N)
        # Check if strain is effectively extinct
        if x[i] <= 1e-10
            # If so x should be set to zero and should not change from that
            dx[i] = 0.0
            x[i] = 0.0
            # In this case the energy concentration should also be fixed to zero
            dx[ps.N + ps.M + i] = 0.0
            x[ps.N + ps.M + i] = 0.0
            # Corresponding proteome fraction also shouldn't shift
            dx[2 * ps.N + ps.M + i] = 0.0
        else
            # find growth rate for strains that aren't extinct
            λ = λs(x[ps.N + ps.M + i], x[2 * ps.N + ps.M + i], ps.mics[i])
            # (growth rate - death rate)*population
            dx[i] = (λ - ps.mics[i].d) * x[i]
            # Now find optimal ribosome fraction
            ϕR = ϕ_R(x[ps.N + ps.M + i], ps.mics[i])
            # This introduces a time delay
            τ = ps.mics[i].fd / λ
            # Then update actual ribosome fraction
            dx[2 * ps.N + ps.M + i] = (ϕR - x[2 * ps.N + ps.M + i]) / τ
            # Energy intake is zero
            J = 0
            # Loop over all reactions to find energy gained by them
            for j in 1:(ps.mics[i].R)
                J += ps.mics[i].η[j] * rate[i, ps.mics[i].Reacs[j]]
            end
            # Add energy intake and subtract translation and dilution from the energy concentration
            dx[ps.N + ps.M + i] = J -
                                  (ps.mics[i].MC * ps.mics[i].ρ + x[ps.N + ps.M + i]) * λ
        end
    end
    # Do basic resource dynamics
    for i in (ps.N + 1):(ps.N + ps.M)
        # fist add external supply of resource and decay
        dx[i] = ps.κ[i - ps.N] - ps.δ[i - ps.N] * x[i]
    end
    # Fine to here
    # Then loop over microbes
    for i in 1:(ps.N)
        # Loop over reactions for specific microbe
        for j in 1:(ps.mics[i].R)
            # Increase the product
            dx[ps.N + ps.reacs[ps.mics[i].Reacs[j]].Prd] += rate[i, ps.mics[i].Reacs[j]] *
                                                            x[i] / NA
            # and decrease the reactant
            dx[ps.N + ps.reacs[ps.mics[i].Reacs[j]].Rct] -= rate[i, ps.mics[i].Reacs[j]] *
                                                            x[i] / NA
        end
    end
    # Final step to correct for any concentrations that have dropped below threshold (1e-15)
    for i in (ps.N + 1):(ps.N + ps.M)
        # If the rate of change is above a threshold (1e-20) they are not altered
        if x[i] < 1e-15
            x[i] = 1e-15
            dx[i] = 0.0
        end
    end
    # Any ATP numbers that have gone below 0.33 should be removed
    for i in (ps.N + ps.M + 1):(2 * ps.N + ps.M)
        if x[i] < 0.33
            x[i] = 0.0
            dx[i] = 0.0
        end
    end
    return (dx)
end

# Simulation code to run one instance of the simulation with a user defined starting condition
# ps is parameter set, Tmax is the time to integrate to, pop, conc, as and ϕs are the initial conditions
function full_simulate(ps::FullParameters, Tmax::Float64, pop::Array{Float64, 1},
                       conc::Array{Float64, 1},
                       as::Array{Float64, 1}, ϕs::Array{Float64, 1})
    @assert length(pop)==ps.N "From parameter set expected $(ps.N) strains"
    @assert length(conc)==ps.M "From parameter set expected $(ps.M) metabolites"
    @assert length(as)==length(pop) "Every strain must have an energy concentration"
    @assert length(ϕs)==length(pop) "Every strain must have a ribosome fraction"
    # Preallocate memory
    rate = zeros(ps.N, ps.O)
    # Now substitute preallocated memory in
    dyns!(dx, x, ps, t) = full_dynamics!(dx, x, ps, rate, t)
    # Find time span for this step
    tspan = (0, Tmax)
    x0 = [pop; conc; as; ϕs]
    # Then setup and solve the problem
    prob = ODEProblem(dyns!, x0, tspan, ps)
    # Still generates problems, not sure if I have to change a solver option or what
    sol = DifferentialEquations.solve(prob)
    return (sol', sol.t)
end

# Same as the above but with a specific strain given a massively increased death rate
# This should help us investigate syntrophic pairs
function full_simulate_syn(ps::FullParameters, Tmax::Float64, pop::Array{Float64, 1},
                           conc::Array{Float64, 1},
                           as::Array{Float64, 1}, ϕs::Array{Float64, 1},
                           inds::Array{Int64, 1})
    @assert length(pop)==ps.N "From parameter set expected $(ps.N) strains"
    @assert length(conc)==ps.M "From parameter set expected $(ps.M) metabolites"
    @assert length(as)==length(pop) "Every strain must have an energy concentration"
    @assert length(ϕs)==length(pop) "Every strain must have a ribosome fraction"
    # Preallocate memory
    rate = zeros(ps.N, ps.O)
    # Extract old death rate
    d = ps.mics[1].d
    # Increase by a factor of 10
    nd = 25.0 * d
    # Preallocate new microbes
    nmics = Array{MicrobeP, 1}(undef, length(ps.mics))
    # Loop over old microbes
    for i in 1:length(ps.mics)
        # Check if any change needs to be made to the microbes
        if i ∉ inds
            nmics[i] = ps.mics[i]
        else
            # Extract old microbe
            m = ps.mics[i]
            # Make new microbe
            nmics[i] = make_MicrobeP(m.MC, m.γm, m.ρ, m.Kγ, m.Pb, nd, m.ϕH, m.KΩ, m.fd, m.R,
                                     m.Reacs, m.η, m.kc, m.KS, m.kr, m.n, m.ϕP)
        end
    end
    # Use to make a new parameter set
    ps2 = make_FullParameters(ps.N, ps.M, ps.O, ps.T, ps.κ, ps.δ, ps.reacs, nmics)
    # Now substitute preallocated memory in
    dyns!(dx, x, ps2, t) = full_dynamics!(dx, x, ps2, rate, t)
    # Find time span for this step
    tspan = (0, Tmax)
    x0 = [pop; conc; as; ϕs]
    # Then setup and solve the problem
    prob = ODEProblem(dyns!, x0, tspan, ps2)
    # Still generates problems, not sure if I have to change a solver option or what
    sol = DifferentialEquations.solve(prob)
    return (sol', sol.t)
end

# Same dynamics function altered to do detailed testing
function test_dynamics!(dx::Array{Float64, 1}, x::Array{Float64, 1}, ps::FullParameters,
                        rate::Array{Float64, 2}, t::Float64)
    # loop over the reactions to find reaction rate for each reaction for each strain
    for j in 1:(ps.O)
        # Find substrate and product for this reaction
        for i in 1:(ps.N)
            # Check if microbe i performs reaction j
            if j ∈ ps.mics[i].Reacs
                # Find index of this reaction in microbe
                k = findfirst(x -> x == j, ps.mics[i].Reacs)
                # Find amount of enzyme E
                E = Eα(x[2 * ps.N + ps.M + i], ps.mics[i], k)
                # Then finally calculate reaction rate
                rate[i, j] = qs(x[ps.N + ps.reacs[j].Rct], x[ps.N + ps.reacs[j].Prd], E, k,
                                ps.mics[i], ps.T, ps.reacs[ps.mics[i].Reacs[k]])
            else
                rate[i, j] = 0.0
            end
        end
    end
    # Now want to use the rate matrix in the consumer dynamics
    for i in 1:(ps.N)
        # Check if strain is effectively extinct
        if x[i] <= 1e-10
            # If so x should be set to zero and should not change from that
            dx[i] = 0.0
            x[i] = 0.0
            # In this case the energy concentration should also be fixed to zero
            dx[ps.N + ps.M + i] = 0.0
            x[ps.N + ps.M + i] = 0.0
            # Corresponding proteome fraction also shouldn't shift
            dx[2 * ps.N + ps.M + i] = 0.0
        else
            # find growth rate for strains that aren't extinct
            λ = λs(x[ps.N + ps.M + i], x[2 * ps.N + ps.M + i], ps.mics[i])
            # (growth rate - death rate)*population
            dx[i] = (λ - ps.mics[i].d) * x[i]
            # Now find optimal ribosome fraction
            ϕR = ϕ_R(x[ps.N + ps.M + i], ps.mics[i])
            # This introduces a time delay
            τ = ps.mics[i].fd / λ
            # Then update actual ribosome fraction
            dx[2 * ps.N + ps.M + i] = (ϕR - x[2 * ps.N + ps.M + i]) / τ
            # Energy intake is zero
            J = 0
            # Loop over all reactions to find energy gained by them
            for j in 1:(ps.mics[i].R)
                J += ps.mics[i].η[j] * rate[i, ps.mics[i].Reacs[j]]
            end
            # Add energy intake and subtract translation and dilution from the energy concentration
            dx[ps.N + ps.M + i] = J -
                                  (ps.mics[i].MC * ps.mics[i].ρ + x[ps.N + ps.M + i]) * λ
        end
    end
    # Do basic resource dynamics
    for i in (ps.N + 1):(ps.N + ps.M)
        # fist add external supply of resource and decay
        dx[i] = ps.κ[i - ps.N] - ps.δ[i - ps.N] * x[i]
    end
    # Fine to here
    # Then loop over microbes
    for i in 1:(ps.N)
        # Loop over reactions for specific microbe
        for j in 1:(ps.mics[i].R)
            # Increase the product
            dx[ps.N + ps.reacs[ps.mics[i].Reacs[j]].Prd] += rate[i, ps.mics[i].Reacs[j]] *
                                                            x[i] / NA
            # and decrease the reactant
            dx[ps.N + ps.reacs[ps.mics[i].Reacs[j]].Rct] -= rate[i, ps.mics[i].Reacs[j]] *
                                                            x[i] / NA
        end
    end
    # Final step to correct for any concentrations that have dropped below threshold (1e-15)
    for i in (ps.N + 1):(ps.N + ps.M)
        # If the rate of change is above a threshold (1e-20) they are not altered
        if x[i] < 1e-15 && dx[i] <= 1e-20
            x[i] = 0.0
            dx[i] = 0.0
        end
    end
    # Any ATP numbers that have gone below 0.33 should be removed
    for i in (ps.N + ps.M + 1):(2 * ps.N + ps.M)
        if x[i] < 0.33
            x[i] = 0.0
            dx[i] = 0.0
        end
    end
    return (dx)
end

# Same function as above but for detailed testing
function test_full_simulate(ps::FullParameters, Tmax::Float64, pop::Array{Float64, 1},
                            conc::Array{Float64, 1},
                            as::Array{Float64, 1}, ϕs::Array{Float64, 1})
    @assert length(pop)==ps.N "From parameter set expected $(ps.N) strains"
    @assert length(conc)==ps.M "From parameter set expected $(ps.M) metabolites"
    @assert length(as)==length(pop) "Every strain must have an energy concentration"
    @assert length(ϕs)==length(pop) "Every strain must have a ribosome fraction"
    # Preallocate memory
    rate = zeros(ps.N, ps.O)
    # Now substitute preallocated memory in
    dyns!(dx, x, ps, t) = test_dynamics!(dx, x, ps, rate, t)
    # Find time span for this step
    tspan = (0, Tmax)
    x0 = [pop; conc; as; ϕs]
    # Then setup and solve the problem
    println("Test simulation started.")
    prob = ODEProblem(dyns!, x0, tspan, ps)
    sol = DifferentialEquations.solve(prob)
    println(sol.destats) # Also useful
    println(sol.retcode) # Useful
    return (sol', sol.t)
end
