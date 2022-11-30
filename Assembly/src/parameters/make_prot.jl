# Script that makes the parameters needed for simulation of the single species proteome model
# Probably should get merged eventually

export initialise_prot, initialise_prot_M, initialise_prot_T, initialise_prot_fix,
       initialise_prot_KΩ, initialise_prot_gl
export initialise_ρ

# function to generate parameter set for the model with inhibition
function initialise_prot(inhib::Bool)
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
    # Now make the reaction
    ΔG = -6e5 # Relatively small Gibbs free energy change
    r = make_Reaction(1, 1, 2, ΔG)
    # η chosen so that a substantial portion of the Gibbs free energy is retained
    if inhib == false
        η = 0.9 * (-ΔG / ΔGATP)
    else
        η = 1.0 * (-ΔG / ΔGATP)
    end
    # The reversibility factor remains the same as previously
    kr = 10.0
    # From Posfai et al (2017) dilution rate 0.21 per hour
    δ = 6.0e-5 * ones(2) # Metabolite dilution rate
    # Set death rate as equal to the dilution rate, just to control number of parameters really
    d = 6.0e-5
    # Human blood glucose is approximately 5.5 m mol per litre (wikipedia)
    # Sensible order of magnitude to aim for set κ/δ = 5.5e-3
    κ = [3.3e-7, 0.0] # Metabolite supply rate
    # Assume that half saturation occurs at a quarter κ/δ
    KS = (1 / 4) * 5.5e-3
    # From wikipedia an average enzyme has a k to KS ratio of 10^5 M^-1 s^-1
    # This would give us a k of 137.5, sensible to assume an above average rate
    # Though should be reduced by the fact we include uptake as well as metabolism
    # Choosing k = 500 means we match the maximum glucose uptake rate seen in Natarajan et al (2000)
    # of 3*10^7 molecules per second.
    # The above is a sensible argument but 1.0 gives a more reasonable ATP concentration.
    kc = 1.0
    # The number of ATP per translation step, including the cost of amino acid synthesis
    # This figure is taken from Lynch and Marinov 2015
    ρ = 29.0
    # This is currently a parameter which I am fiddling
    Kγ = 5e8
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # Give omega a fairly arbitrary value for now, would be expected to be of similar order to Kγ
    KΩ = 2 * Kγ
    # Number of doublings required to dilute to 1%
    fd = log(100) / log(2)
    # Now make the parameter set
    ps = make_ProtParameters(MC, γm, T, η, KS, kr, kc, ρ, Kγ, d, Pb, ϕH, KΩ, fd, r, n, δ, κ)
    return (ps)
end

# function to generate parameter set for the model with inhibition option to change cell mass
function initialise_prot_M(MC::Int64)
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # Elongation rate also taken from Bremer (1996), convert from minutes to seconds
    γm = 1260.0 / 60.0
    # Now preallocate protein masses
    n = zeros(Int64, 3)
    # ribosome mass taken from Keseler IM, et al. (2011)
    n[1] = 7459
    # Other protein mass averaged from Brandt F, et al. (2009)
    n[2:3] .= 300
    # Now make the reaction
    ΔG = -6e5 # Relatively small Gibbs free energy change
    r = make_Reaction(1, 1, 2, ΔG)
    # η chosen so that a substantial portion of the Gibbs free energy is retained
    η = 0.9 * (-ΔG / ΔGATP)
    # The reversibility factor remains the same as previously
    kr = 10.0
    # From Posfai et al (2017) dilution rate 0.21 per hour
    δ = 6.0e-5 * ones(2) # Metabolite dilution rate
    # Set death rate as equal to the dilution rate, just to control number of parameters really
    d = 6.0e-5
    # Human blood glucose is approximately 5.5 m mol per litre (wikipedia)
    # Sensible order of magnitude to aim for set κ/δ = 5.5e-3
    κ = [3.3e-7, 0.0] # Metabolite supply rate
    # Assume that half saturation occurs at a quarter κ/δ
    KS = (1 / 4) * 5.5e-3
    # From wikipedia an average enzyme has a k to KS ratio of 10^5 M^-1 s^-1
    # This would give us a k of 137.5, sensible to assume an above average rate
    # Though should be reduced by the fact we include uptake as well as metabolism
    # Choosing k = 500 means we match the maximum glucose uptake rate seen in Natarajan et al (2000)
    # of 3*10^7 molecules per second.
    # The above is a sensible argument but 1.0 gives a more reasonable ATP concentration.
    kc = 1.0
    # The number of ATP per translation step, including the cost of amino acid synthesis
    # This figure is taken from Lynch and Marinov 2015
    ρ = 29.0
    # This is currently a parameter which I am fiddling
    Kγ = 5e8 * (MC / (10^8))
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # Give omega a fairly arbitrary value for now, would be expected to be of similar order to Kγ
    KΩ = 2 * Kγ
    # Number of doublings required to dilute to 1%
    fd = log(100) / log(2)
    # Now make the parameter set
    ps = make_ProtParameters(MC, γm, T, η, KS, kr, kc, ρ, Kγ, d, Pb, ϕH, KΩ, fd, r, n, δ, κ)
    return (ps)
end

# function to generate parameter set for the model with inhibition option to change temperature
function initialise_prot_T(T::Float64, inhib::Bool)
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
    # Now make the reaction
    ΔG = -6e5 # Relatively small Gibbs free energy change
    r = make_Reaction(1, 1, 2, ΔG)
    # η chosen so that a substantial portion of the Gibbs free energy is retained
    # η chosen so that a substantial portion of the Gibbs free energy is retained
    if inhib == false
        η = 0.9 * (-ΔG / ΔGATP)
    else
        η = 1.1 * (-ΔG / ΔGATP)
    end
    # The reversibility factor remains the same as previously
    kr = 10.0
    # From Posfai et al (2017) dilution rate 0.21 per hour
    δ = 6.0e-5 * ones(2) # Metabolite dilution rate
    # Set death rate as equal to the dilution rate, just to control number of parameters really
    d = 6.0e-5
    # Human blood glucose is approximately 5.5 m mol per litre (wikipedia)
    # Sensible order of magnitude to aim for set κ/δ = 5.5e-3
    κ = [3.3e-7, 0.0] # Metabolite supply rate
    # Assume that half saturation occurs at a quarter κ/δ
    KS = (1 / 4) * 5.5e-3
    # From wikipedia an average enzyme has a k to KS ratio of 10^5 M^-1 s^-1
    # This would give us a k of 137.5, sensible to assume an above average rate
    # Though should be reduced by the fact we include uptake as well as metabolism
    # Choosing k = 500 means we match the maximum glucose uptake rate seen in Natarajan et al (2000)
    # of 3*10^7 molecules per second.
    # The above is a sensible argument but 1.0 gives a more reasonable ATP concentration.
    kc = 1.0
    # The number of ATP per translation step, including the cost of amino acid synthesis
    # This figure is taken from Lynch and Marinov 2015
    ρ = 29.0
    # This is currently a parameter which I am fiddling
    Kγ = 5e8 * (MC / (10^8))
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # Give omega a fairly arbitrary value for now, would be expected to be of similar order to Kγ
    KΩ = 2 * Kγ
    # Number of doublings required to dilute to 1%
    fd = log(100) / log(2)
    # Now make the parameter set
    ps = make_ProtParameters(MC, γm, T, η, KS, kr, kc, ρ, Kγ, d, Pb, ϕH, KΩ, fd, r, n, δ, κ)
    return (ps)
end

# function to generate parameter set for the model outside the chemostat
function initialise_prot_fix(Kγ::Float64, KΩ::Float64, kc::Float64, η::Float64 = 7.2)
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
    # Now make the reaction
    ΔG = -6e5 # Relatively small Gibbs free energy change
    r = make_Reaction(1, 1, 2, ΔG)
    # # η chosen so that a substantial portion of the Gibbs free energy is retained
    # η = 0.9*(-ΔG/ΔGATP)
    # The reversibility factor remains the same as previously
    kr = 10.0
    # In this case there is (effectively) no removal
    δ = 1e-9 * ones(2)
    # Set death rate to (basically) zero
    d = 6.0e-5
    # Also no supply in this case
    κ = [0.0, 0.0] # Metabolite supply rate
    # Assume that half saturation occurs at a quarter κ/δ
    KS = (1 / 4) * 5.5e-3
    # The number of ATP per translation step, including the cost of amino acid synthesis
    # This figure is taken from Lynch and Marinov 2015
    ρ = 29.0
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # Number of doublings required to dilute to 1%
    fd = log(100) / log(2)
    # Now make the parameter set
    ps = make_ProtParameters(MC, γm, T, η, KS, kr, kc, ρ, Kγ, d, Pb, ϕH, KΩ, fd, r, n, δ, κ)
    return (ps)
end

# function to generate parameter set for the model with inhibition
function initialise_prot_KΩ(ps::ProtParameters, KΩ::Float64)
    # Now make the parameter set
    ps = make_ProtParameters(ps.MC, ps.γm, ps.T, ps.η, ps.KS, ps.kr, ps.kc, ps.ρ, ps.Kγ,
                             ps.d, ps.Pb, ps.ϕH, KΩ, ps.fd, ps.r, ps.n, ps.δ, ps.κ)
    return (ps)
end

# function to generate parameter set for the model with inhibition
function initialise_prot_gl(γm::Float64, ΔG::Float64)
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # Cell mass is taken from Bremer H, Dennis P (1996) Modulation of chemical
    # composition and other parameters of the cell by growth rate (Book chapter).
    MC = 10^8
    # Now preallocate protein masses
    n = zeros(Int64, 3)
    # ribosome mass taken from Keseler IM, et al. (2011)
    n[1] = 7459
    # Other protein mass averaged from Brandt F, et al. (2009)
    n[2:3] .= 300
    # Now make the reaction
    r = make_Reaction(1, 1, 2, ΔG)
    # η chosen so that a substantial portion of the Gibbs free energy is retained
    η = 0.9 * (-ΔG / ΔGATP)
    # The reversibility factor remains the same as previously
    kr = 10.0
    # Considering batch culture so want this to be basically zero
    δ = 1.0e-99 * ones(2) # Metabolite dilution rate
    # Set death rate as equal to the dilution rate, just to control number of parameters really
    d = 6.0e-5
    # Human blood glucose is approximately 5.5 m mol per litre (wikipedia)
    # Sensible order of magnitude to aim for set κ/δ = 5.5e-3
    κ = [3.3e-7, 0.0] # Metabolite supply rate
    # Assume that half saturation occurs at a quarter κ/δ
    KS = (1 / 4) * 5.5e-3
    # From wikipedia an average enzyme has a k to KS ratio of 10^5 M^-1 s^-1
    # This would give us a k of 137.5, sensible to assume an above average rate
    # Though should be reduced by the fact we include uptake as well as metabolism
    # Choosing k = 500 means we match the maximum glucose uptake rate seen in Natarajan et al (2000)
    # of 3*10^7 molecules per second.
    # The above is a sensible argument but 1.0 gives a more reasonable ATP concentration.
    kc = 1.0
    # The number of ATP per translation step, including the cost of amino acid synthesis
    # This figure is taken from Lynch and Marinov 2015
    ρ = 29.0
    # This is currently a parameter which I am fiddling
    Kγ = 5e8
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # Give omega a fairly arbitrary value for now, would be expected to be of similar order to Kγ
    KΩ = 2 * Kγ
    # Number of doublings required to dilute to 1%
    fd = log(100) / log(2)
    # Now make the parameter set
    ps = make_ProtParameters(MC, γm, T, η, KS, kr, kc, ρ, Kγ, d, Pb, ϕH, KΩ, fd, r, n, δ, κ)
    return (ps)
end

# function to make parameter set for ρ vs γmax tradeoff
function initialise_ρ(fac::Float64)
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
    # Now make the reaction
    ΔG = -6e5 # Relatively small Gibbs free energy change
    r = make_Reaction(1, 1, 2, ΔG)
    # η chosen so that a substantial portion of the Gibbs free energy is retained
    η = 0.9 * (-ΔG / ΔGATP)
    # The reversibility factor remains the same as previously
    kr = 10.0
    # From Posfai et al (2017) dilution rate 0.21 per hour
    δ = 6.0e-5 * ones(2) # Metabolite dilution rate
    # Set death rate as equal to the dilution rate, just to control number of parameters really
    d = 6.0e-5
    # Human blood glucose is approximately 5.5 m mol per litre (wikipedia)
    # Sensible order of magnitude to aim for set κ/δ = 5.5e-3
    κ = [3.3e-7, 0.0] # Metabolite supply rate
    # Assume that half saturation occurs at a quarter κ/δ
    KS = (1 / 4) * 5.5e-3
    # From wikipedia an average enzyme has a k to KS ratio of 10^5 M^-1 s^-1
    # This would give us a k of 137.5, sensible to assume an above average rate
    # Though should be reduced by the fact we include uptake as well as metabolism
    # Choosing k = 500 means we match the maximum glucose uptake rate seen in Natarajan et al (2000)
    # of 3*10^7 molecules per second.
    # The above is a sensible argument but 1.0 gives a more reasonable ATP concentration.
    kc = 1.0
    # The number of ATP per translation step, including the cost of amino acid synthesis
    # This figure is taken from Lynch and Marinov 2015
    ρ = 29.0
    # This is currently a parameter which I am fiddling
    Kγ = 5e8
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # Give omega a fairly arbitrary value for now, would be expected to be of similar order to Kγ
    KΩ = 2 * Kγ
    # Number of doublings required to dilute to 1%
    fd = log(100) / log(2)
    # Need to now to the rescaling
    γm = γm * fac
    ρ = ρ * fac
    # Now make the parameter set
    ps = make_ProtParameters(MC, γm, T, η, KS, kr, kc, ρ, Kγ, d, Pb, ϕH, KΩ, fd, r, n, δ, κ)
    return (ps)
end
