# A script to analyse single population growth to test and parameterise the new model
using TradeOff
using Random

# Simply want to test if the new function I've written works
function test()
    println("Compiled and starting simulation")
    # Two metabolites 1 reaction for the sake of simplicity
    M = 2
    O = 1
    # First generate parameter set
    ps = initialise(M,O)
    # Only one reaction so upper and lower bound have to be one
    Rl = 1
    Ru = 1
    # Choose new parameter values
    μ = 1e-1
    ω = 0.5
    KΩ = 1e9
    Kχ = 1e9
    # Then make a microbes
    mic = new_mic(M,μ,ω,KΩ,Kχ,Rl,Ru)
    println(mic)
    return(nothing)
end

@time test()
