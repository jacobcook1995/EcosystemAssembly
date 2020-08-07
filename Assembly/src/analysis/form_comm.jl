# A script to form stable communities using the full model that then saves them for later use
using Assembly

# function to test that the new stuff I'm writing actually works
function test()
    println("Compiled")
    # Assume that half saturation occurs at a quarter κ/δ
    KS = (1/4)*5.5e-3
    # From wikipedia an average enzyme has a k to KS ratio of 10^5 M^-1 s^-1
    # This would give us a k of 137.5, sensible to assume an above average rate
    # Though should be reduced by the fact we include uptake as well as metabolism
    # Choosing k = 500 means we match the maximum glucose uptake rate seen in Natarajan et al (2000)
    # of 3*10^7 molecules per second.
    # The above is a sensible argument but 1.0 gives a more reasonable ATP concentration.
    kc = 1.0
    # The reversibility factor remains the same as previously
    kr = 10.0
    # Assume microbes have 2 reactions each
    mR = 2.0
    sdR = 0.0
    # Case of 8 metabolites and 1 strain
    N = 1
    M = 8
    # Use formula to find how many reactions this implies
    O = 2*M - 3
    ps = initialise(N,M,O::Int64,mR,sdR,kc,KS,kr)
    println(ps)
    return(nothing)
end

@time test()
