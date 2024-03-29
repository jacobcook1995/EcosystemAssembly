# Function in here are designed to set up chemostat runs to test impact of
# parametrisations on growth rates
export ω_test

# function to implement the consumer resource dynamics
function chemo_dynamics!(dx::Array{Float64, 1},
                         x::Array{Float64, 1},
                         ms::Array{Microbe, 1},
                         ps::TOParameters,
                         C::Array{Float64, 1},
                         rate::Array{Float64, 2},
                         t::Float64)
    # loop over the reactions to find reaction rate for each reaction for each strain
    for j in 1:(ps.O)
        # Find substrate and product for this reaction
        for i in 1:length(ms)
            # Check if microbe i performs reaction j
            if j ∈ ms[i].Reacs
                # Find index of this reaction in microbe
                k = findfirst(x -> x == j, ms[i].Reacs)
                # Find amount of enzyme E
                E = Eα(x[2 * length(ms) + i], ms[i], k)
                # Then finally calculate reaction rate
                rate[i, j] = qs(C[ps.reacs[j].Rct],
                                C[ps.reacs[j].Prd],
                                E,
                                k,
                                ms[i],
                                ps.T,
                                ps.reacs[ms[i].Reacs[k]])
            else
                rate[i, j] = 0.0
            end
        end
    end
    # Now want to use the rate matrix in the consumer dynamics
    for i in 1:length(ms)
        # Check if strain is effectively extinct
        if x[i] <= 1e-5
            # If so x should be set to zero and should not change from that
            dx[i] = 0.0
            x[i] = 0.0
            # In this case the energy concentration should also be fixed to zero
            dx[length(ms) + i] = 0.0
            x[length(ms) + i] = 0.0
            # Corresponding proteome fraction also shouldn't shift
            dx[2 * length(ms) + i] = 0.0
        else
            # find growth rate for strains that aren't extinct
            λ = λs(x[length(ms) + i], x[2 * length(ms) + i], ms[i])
            # (growth rate - death rate)*population
            dx[i] = (λ - ms[i].d) * x[i]
            # Now find optimal ribosome fraction
            ϕR = ϕ_R(x[length(ms) + i], ms[i])
            # This introduces a time delay
            τ = ms[i].fd / λ
            # Then update actual ribosome fraction
            dx[2 * length(ms) + i] = (ϕR - x[2 * length(ms) + i]) / τ
            # Energy intake is zero
            J = 0
            # Loop over all reactions to find energy gained by them
            for j in 1:(ms[i].R)
                J += ms[i].η[j] * rate[i, ms[i].Reacs[j]]
            end
            # Add energy intake and subtract translation and dilution from the energy concentration
            dx[length(ms) + i] = J - (ms[i].MC * ms[i].χl + x[length(ms) + i]) * λ
        end
    end
    # Any ATP numbers that have gone below 0.33 should be removed
    for i in (length(ms) + 1):(2 * length(ms))
        if x[i] < 0.33
            x[i] = 0.0
            dx[i] = 0.0
        end
    end
    return (dx)
end

# function to test for single population growth
function chemo(ps::TOParameters,
               pop::Float64,
               conc::Float64,
               as::Float64,
               ϕs::Float64,
               mic::Microbe,
               Tmax::Float64)
    # Preallocate memory
    rate = zeros(1, ps.O)
    # Set constant concentration
    concs = [conc, 1e-7]
    # Now substitute preallocated memory in
    dyns!(dx, x, ms, t) = chemo_dynamics!(dx, x, ms, ps, concs, rate, t)
    # Find time span for this step
    tspan = (0, Tmax)
    # Make appropriate initial condition
    x0 = [pop; as; ϕs]
    # Then setup and solve the problem
    prob = ODEProblem(dyns!, x0, tspan, [mic])
    sol = DifferentialEquations.solve(prob)
    return (sol', sol.t)
end

function ω_test()
    println("Compiled")
    # Define basic simulation parameters
    M = 2
    d = 6e-5
    μrange = 5e7 # Set high to avoid thermodynamic inhibition
    # Predefine a set of omega values
    ωs = collect(range(0.1, 0.9, length = 19))
    # First define fixed microbial variables
    η = 3.0 # Reasonably high
    PID = randstring(['0':'9'; 'a':'f'])
    MC = 10^8
    nr = 7459
    ns = 300
    γm = 1260.0 / 60.0
    Kγ = 5e8
    χl = 29.0
    Pb = 0.7
    ϕH = 0.45
    fd = log(100) / log(2)
    KS = (1 / 4) * 5.5e-3
    kc = 10.0
    kr = 10.0
    # Use formula to calculate how many reactions are implied
    O = floor(Int64, M * (M - 1) / 2)
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # Make parameter set
    ps = initialise(M, O, μrange)
    # Generate fixed reaction
    RP, ΔG = fix_reactions(O, M, μrange, T)
    # Preallocate vector to store single reaction
    reacs = Array{Reaction, 1}(undef, O)
    for i in 1:O
        reacs[i] = make_Reaction(i, RP[i, 1], RP[i, 2], ΔG[i])
    end
    # Only one reaction, which all species possess
    R = 1
    Reacs = [1]
    # Now preallocate protein masses
    n = zeros(Int64, 2 + R)
    # First element is ribosome mass
    n[1] = nr
    # Second is housekeeping
    n[2] = ns
    # Determine the others based on reactions
    for j in 1:R
        n[2 + j] = ns * (reacs[Reacs[j]].Prd - reacs[Reacs[j]].Rct)
    end
    # Preallocate array of fixed microbes
    fix = Array{Microbe, 1}(undef, length(ωs))
    var = Array{Microbe, 1}(undef, length(ωs))
    # Set fixed Ω
    Ωf = 2.5e8
    for i in eachindex(ωs)
        # Can finally generate microbe
        fix[i] = make_Microbe(MC,
                              γm,
                              Kγ,
                              χl,
                              Pb,
                              d,
                              ϕH,
                              Ωf,
                              fd,
                              ωs[i],
                              R,
                              Reacs,
                              [η],
                              [kc],
                              [KS],
                              [kr],
                              n,
                              [1.0],
                              i,
                              PID)
        var[i] = make_Microbe(MC,
                              γm,
                              Kγ,
                              χl,
                              Pb,
                              d,
                              ϕH,
                              Ωf,
                              fd,
                              ωs[i],
                              R,
                              Reacs,
                              [η],
                              [kc],
                              [KS],
                              [kr],
                              n,
                              [1.0],
                              i + length(ωs),
                              PID)
    end
    # Choose sensible initial values
    ϕi = 0.01 # Start at low value
    ai = 1e5
    Cl = 0.001
    Ch = 1.0
    Ni = 1e-3
    # Preallocate containers
    λlf = zeros(length(ωs))
    λhf = zeros(length(ωs))
    λlv = zeros(length(ωs))
    λhv = zeros(length(ωs))
    # Choose simulation window
    Tmax = 5e5
    # Loop over all the fixed species
    for i in eachindex(ωs)
        # Simulate each population (first for low conc)
        C, T = chemo(ps, Ni, Cl, ai, ϕi, fix[i], Tmax)
        # Find growth rate of final point
        λlf[i] = λs(C[end, 2], C[end, 3], fix[i])
        # Then simulate the high conc case
        C, T = chemo(ps, Ni, Ch, ai, ϕi, fix[i], Tmax)
        # Find growth rate of final point
        λhf[i] = λs(C[end, 2], C[end, 3], fix[i])
        # Then move onto the variable cases
        C, T = chemo(ps, Ni, Cl, ai, ϕi, var[i], Tmax)
        # Find growth rate of final point
        λlv[i] = λs(C[end, 2], C[end, 3], var[i])
        # Then simulate the high conc case
        C, T = chemo(ps, Ni, Ch, ai, ϕi, var[i], Tmax)
        # Find growth rate of final point
        λhv[i] = λs(C[end, 2], C[end, 3], var[i])
    end
    return (ωs, λlf, λhf, λlv, λhv)
end
