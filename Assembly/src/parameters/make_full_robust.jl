# Script that makes the parameters needed for simulation of the full proteome model
export options_titles, initialise_np, initialise_Pb, initialise_ϕH, initialise_sat

# Function to provide names for each robustness option
function options_titles(opt::Int64)
    if opt == 1
        return ("ChangeMetMass")
    elseif opt == 2
        return ("ChangeBinding")
    elseif opt == 3
        return ("ChangeHouse")
    elseif opt == 4
        return ("RaiseGamma")
    elseif opt == 5
        return ("HighSat")
    else
        return ("LowSat")
    end
end

# function to generate parameter set for the model with inhibition
function initialise_np(N::Int64, M::Int64, O::Int64, Rl::Int64, Ru::Int64, kc::Float64,
                       KS::Float64,
                       kr::Float64, syn::Bool, μrange::Float64, np::Int64)
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # Cell mass is taken from Bremer H, Dennis P (1996) Modulation of chemical
    # composition and other parameters of the cell by growth rate (Book chapter).
    MC = 10^8
    # Elongation rate also taken from Bremer (1996), convert from minutes to seconds
    γm = 1260.0 / 60.0
    # Now preallocate protein masses
    n = zeros(Int64, 3)
    # ribosome mass taken from Keseler IM, et al. (2011)
    n[1] = 7459
    # In this case the metabolic protein mass is allowed to vary
    n[2:3] .= np
    # Need to think carefully about what a reasonable parameter value is here
    # This has a large impact on the fraction, so should be careful with this one
    d = 6.0e-5
    # The number of ATP per translation step, including the cost of amino acid synthesis
    # This figure is taken from Lynch and Marinov 2015
    ρ = 29.0
    # Back to old arbitrary figures, think this might be the best route
    Kγ = 5e8
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # Back to old arbitrary figures, think this might be the best route
    KΩ = 1e9
    # Number of doublings required to dilute to 1%
    fd = log(100) / log(2)
    # From Posfai et al (2017) dilution rate 0.21 per hour
    δ = 6.0e-5 * ones(M) # Metabolite dilution rate
    # Human blood glucose is approximately 5.5 m mol per litre (wikipedia)
    # Sensible order of magnitude to aim for set κ/δ = 5.5e-3
    κ = zeros(M)
    # All but resource 1 is not supplied
    κ[1] = 3.3e-7 # Metabolite supply rate
    # Generate fixed set of reactions
    RP, ΔG = fix_reactions(O, M, μrange, T)
    # Preallocate vector of reactions
    reacs = Array{Reaction, 1}(undef, O)
    for i in 1:O
        reacs[i] = make_Reaction(i, RP[i, 1], RP[i, 2], ΔG[i])
    end
    # Preallocate vector of microbes
    mics = Array{MicrobeP, 1}(undef, N)
    # Then construct microbes
    for i in 1:N
        # For each microbe generate random set of reactions
        R, Reacs = choose_reactions(O, Rl, Ru)
        # Make vectors of the (fixed) kinetic parameters
        kcs, KSs, krs = kin_rand(kc, KS, kr, R)
        # Reactions given random proportional weightings, done this in the simplest way possible
        ϕP = rand(R)
        ϕP = ϕP / sum(ϕP)
        # Find corresponding η's for these reactions
        η = choose_η_mix(reacs, Reacs, T, syn)
        # Can finally generate microbe
        mics[i] = make_MicrobeP(MC, γm, ρ, Kγ, Pb, d, ϕH, KΩ, fd, R, Reacs, η, kcs, KSs,
                                krs, n, ϕP)
    end
    # Now make the parameter set
    ps = make_FullParameters(N, M, O, T, κ, δ, reacs, mics)
    return (ps)
end

# function to generate parameter set for the model with inhibition
function initialise_Pb(N::Int64, M::Int64, O::Int64, Rl::Int64, Ru::Int64, kc::Float64,
                       KS::Float64,
                       kr::Float64, syn::Bool, μrange::Float64, Pb::Float64)
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # Cell mass is taken from Bremer H, Dennis P (1996) Modulation of chemical
    # composition and other parameters of the cell by growth rate (Book chapter).
    MC = 10^8
    # Elongation rate also taken from Bremer (1996), convert from minutes to seconds
    γm = 1260.0 / 60.0
    # Now preallocate protein masses
    n = zeros(Int64, 3)
    # ribosome mass taken from Keseler IM, et al. (2011)
    n[1] = 7459
    # Other protein mass averaged from Brandt F, et al. (2009)
    n[2:3] .= 300
    # Need to think carefully about what a reasonable parameter value is here
    # This has a large impact on the fraction, so should be careful with this one
    d = 6.0e-5
    # The number of ATP per translation step, including the cost of amino acid synthesis
    # This figure is taken from Lynch and Marinov 2015
    ρ = 29.0
    # Back to old arbitrary figures, think this might be the best route
    Kγ = 5e8
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # Back to old arbitrary figures, think this might be the best route
    KΩ = 1e9
    # Number of doublings required to dilute to 1%
    fd = log(100) / log(2)
    # From Posfai et al (2017) dilution rate 0.21 per hour
    δ = 6.0e-5 * ones(M) # Metabolite dilution rate
    # Human blood glucose is approximately 5.5 m mol per litre (wikipedia)
    # Sensible order of magnitude to aim for set κ/δ = 5.5e-3
    κ = zeros(M)
    # All but resource 1 is not supplied
    κ[1] = 3.3e-7 # Metabolite supply rate
    # Generate fixed set of reactions
    RP, ΔG = fix_reactions(O, M, μrange, T)
    # Preallocate vector of reactions
    reacs = Array{Reaction, 1}(undef, O)
    for i in 1:O
        reacs[i] = make_Reaction(i, RP[i, 1], RP[i, 2], ΔG[i])
    end
    # Preallocate vector of microbes
    mics = Array{MicrobeP, 1}(undef, N)
    # Then construct microbes
    for i in 1:N
        # For each microbe generate random set of reactions
        R, Reacs = choose_reactions(O, Rl, Ru)
        # Make vectors of the (fixed) kinetic parameters
        kcs, KSs, krs = kin_rand(kc, KS, kr, R)
        # Reactions given random proportional weightings, done this in the simplest way possible
        ϕP = rand(R)
        ϕP = ϕP / sum(ϕP)
        # Find corresponding η's for these reactions
        η = choose_η_mix(reacs, Reacs, T, syn)
        # Can finally generate microbe
        mics[i] = make_MicrobeP(MC, γm, ρ, Kγ, Pb, d, ϕH, KΩ, fd, R, Reacs, η, kcs, KSs,
                                krs, n, ϕP)
    end
    # Now make the parameter set
    ps = make_FullParameters(N, M, O, T, κ, δ, reacs, mics)
    return (ps)
end

# function to generate parameter set for the model with inhibition
function initialise_ϕH(N::Int64, M::Int64, O::Int64, Rl::Int64, Ru::Int64, kc::Float64,
                       KS::Float64,
                       kr::Float64, syn::Bool, μrange::Float64, ϕH::Float64)
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # Cell mass is taken from Bremer H, Dennis P (1996) Modulation of chemical
    # composition and other parameters of the cell by growth rate (Book chapter).
    MC = 10^8
    # Elongation rate also taken from Bremer (1996), convert from minutes to seconds
    γm = 1260.0 / 60.0
    # Now preallocate protein masses
    n = zeros(Int64, 3)
    # ribosome mass taken from Keseler IM, et al. (2011)
    n[1] = 7459
    # Other protein mass averaged from Brandt F, et al. (2009)
    n[2:3] .= 300
    # Need to think carefully about what a reasonable parameter value is here
    # This has a large impact on the fraction, so should be careful with this one
    d = 6.0e-5
    # The number of ATP per translation step, including the cost of amino acid synthesis
    # This figure is taken from Lynch and Marinov 2015
    ρ = 29.0
    # Back to old arbitrary figures, think this might be the best route
    Kγ = 5e8
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Back to old arbitrary figures, think this might be the best route
    KΩ = 1e9
    # Number of doublings required to dilute to 1%
    fd = log(100) / log(2)
    # From Posfai et al (2017) dilution rate 0.21 per hour
    δ = 6.0e-5 * ones(M) # Metabolite dilution rate
    # Human blood glucose is approximately 5.5 m mol per litre (wikipedia)
    # Sensible order of magnitude to aim for set κ/δ = 5.5e-3
    κ = zeros(M)
    # All but resource 1 is not supplied
    κ[1] = 3.3e-7 # Metabolite supply rate
    # Generate fixed set of reactions
    RP, ΔG = fix_reactions(O, M, μrange, T)
    # Preallocate vector of reactions
    reacs = Array{Reaction, 1}(undef, O)
    for i in 1:O
        reacs[i] = make_Reaction(i, RP[i, 1], RP[i, 2], ΔG[i])
    end
    # Preallocate vector of microbes
    mics = Array{MicrobeP, 1}(undef, N)
    # Then construct microbes
    for i in 1:N
        # For each microbe generate random set of reactions
        R, Reacs = choose_reactions(O, Rl, Ru)
        # Make vectors of the (fixed) kinetic parameters
        kcs, KSs, krs = kin_rand(kc, KS, kr, R)
        # Reactions given random proportional weightings, done this in the simplest way possible
        ϕP = rand(R)
        ϕP = ϕP / sum(ϕP)
        # Find corresponding η's for these reactions
        η = choose_η_mix(reacs, Reacs, T, syn)
        # Can finally generate microbe
        mics[i] = make_MicrobeP(MC, γm, ρ, Kγ, Pb, d, ϕH, KΩ, fd, R, Reacs, η, kcs, KSs,
                                krs, n, ϕP)
    end
    # Now make the parameter set
    ps = make_FullParameters(N, M, O, T, κ, δ, reacs, mics)
    return (ps)
end

# function to generate parameter set for the model with inhibition
function initialise_sat(N::Int64, M::Int64, O::Int64, Rl::Int64, Ru::Int64, kc::Float64,
                        KS::Float64,
                        kr::Float64, syn::Bool, μrange::Float64, Kγ::Float64, KΩ::Float64)
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # Cell mass is taken from Bremer H, Dennis P (1996) Modulation of chemical
    # composition and other parameters of the cell by growth rate (Book chapter).
    MC = 10^8
    # Elongation rate also taken from Bremer (1996), convert from minutes to seconds
    γm = 1260.0 / 60.0
    # Now preallocate protein masses
    n = zeros(Int64, 3)
    # ribosome mass taken from Keseler IM, et al. (2011)
    n[1] = 7459
    # Other protein mass averaged from Brandt F, et al. (2009)
    n[2:3] .= 300
    # Need to think carefully about what a reasonable parameter value is here
    # This has a large impact on the fraction, so should be careful with this one
    d = 6.0e-5
    # The number of ATP per translation step, including the cost of amino acid synthesis
    # This figure is taken from Lynch and Marinov 2015
    ρ = 29.0
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # Number of doublings required to dilute to 1%
    fd = log(100) / log(2)
    # From Posfai et al (2017) dilution rate 0.21 per hour
    δ = 6.0e-5 * ones(M) # Metabolite dilution rate
    # Human blood glucose is approximately 5.5 m mol per litre (wikipedia)
    # Sensible order of magnitude to aim for set κ/δ = 5.5e-3
    κ = zeros(M)
    # All but resource 1 is not supplied
    κ[1] = 3.3e-7 # Metabolite supply rate
    # Generate fixed set of reactions
    RP, ΔG = fix_reactions(O, M, μrange, T)
    # Preallocate vector of reactions
    reacs = Array{Reaction, 1}(undef, O)
    for i in 1:O
        reacs[i] = make_Reaction(i, RP[i, 1], RP[i, 2], ΔG[i])
    end
    # Preallocate vector of microbes
    mics = Array{MicrobeP, 1}(undef, N)
    # Then construct microbes
    for i in 1:N
        # For each microbe generate random set of reactions
        R, Reacs = choose_reactions(O, Rl, Ru)
        # Make vectors of the (fixed) kinetic parameters
        kcs, KSs, krs = kin_rand(kc, KS, kr, R)
        # Reactions given random proportional weightings, done this in the simplest way possible
        ϕP = rand(R)
        ϕP = ϕP / sum(ϕP)
        # Find corresponding η's for these reactions
        η = choose_η_mix(reacs, Reacs, T, syn)
        # Can finally generate microbe
        mics[i] = make_MicrobeP(MC, γm, ρ, Kγ, Pb, d, ϕH, KΩ, fd, R, Reacs, η, kcs, KSs,
                                krs, n, ϕP)
    end
    # Now make the parameter set
    ps = make_FullParameters(N, M, O, T, κ, δ, reacs, mics)
    return (ps)
end
