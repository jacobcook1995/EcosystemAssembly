# A script to form stable communities using the full model that then saves them for later use
using Assembly

# function to test that the new stuff I'm writing actually works
function test()
    println("Compiled")
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # Cell mass is taken from Bremer H, Dennis P (1996) Modulation of chemical
    # composition and other parameters of the cell by growth rate (Book chapter).
    MC = 10^8
    # Elongation rate also taken from Bremer (1996), convert from minutes to seconds
    γm = 1260.0/60.0
    # Now preallocate protein masses
    n = zeros(Int64,3)
    # ribosome mass taken from Keseler IM, et al. (2011)
    n[1] = 7459
    # Other protein mass averaged from Brandt F, et al. (2009)
    n[2:3] .= 300
    # Set death rate as equal to the dilution rate, just to control number of parameters really
    d = 6.0e-5
    # The number of ATP per translation step, including the cost of amino acid sythesis
    # This figure is taken from Lynch and Marinov 2015
    ρ = 29.0
    # This is currently a paramter which I am fiddling
    Kγ = 5e8
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # Give omega a fairly arbitary value for now, would be expected to be of similar order to Kγ
    KΩ = 2*Kγ
    # Number of doublings required to dilute to 1%
    fd = log(100)/log(2)
    # For simple case only want to consider 1 reaction
    R = 1
    # Now make the reaction
    ΔG = -6e5 # Relatively small Gibbs free energy change
    r = make_Reaction(1,1,2,ΔG)
    # Microbe uses first (and only) reaction
    Reacs = [1]
    # Now make the vector quantities
    # η chosen so that a substantial portion of the Gibbs free energy is retained
    η = [0.9*(-ΔG/ΔGATP)]
    # Assume that half saturation occurs at a quarter κ/δ
    KS = [(1/4)*5.5e-3]
    # From wikipedia an average enzyme has a k to KS ratio of 10^5 M^-1 s^-1
    # This would give us a k of 137.5, sensible to assume an above average rate
    # Though should be reduced by the fact we include uptake as well as metabolism
    # Choosing k = 500 means we match the maximum glucose uptake rate seen in Natarajan et al (2000)
    # of 3*10^7 molecules per second.
    # The above is a sensible argument but 1.0 gives a more reasonable ATP concentration.
    kc = [1.0]
    # The reversibility factor remains the same as previously
    kr = [10.0]
    # First make a microbe
    mic = make_MicrobeP(MC,γm,ρ,Kγ,Pb,d,ϕH,KΩ,fd,R,Reacs,η,kc,KS,kr,n)
    # vectorise this microbe
    mics = [mic]
    # vectorise reaction
    reacs = [r]
    # 1 microbe
    N = 1
    # 2 metabolites
    M = 2
    # 1 reaction
    O = 1
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # From Posfai et al (2017) dilution rate 0.21 per hour
    δ = 6.0e-5*ones(2) # Metabolite dilution rate
    # Human blood glucose is approximatly 5.5 m mol per litre (wikipedia)
    # Sensible order of magnitude to aim for set κ/δ = 5.5e-3
    κ = [3.3e-7,0.0] # Metabolite supply rate
    ps = make_FullParameters(N,M,O,T,κ,δ,reacs,mics)
    return(nothing)
end

@time test()
