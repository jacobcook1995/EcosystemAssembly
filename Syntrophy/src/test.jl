using Syntrophy
# This is a script to write my testing code into
# Anything reusable should be moved into Syntrophy.jl as a seperate callable function

function main()
    # Need to find a suitable testing example tomorrow
    # T = 312.0 # Standard temp
    # sub = [ 2.00*10^(-9), 5.9*10^(-5)]
    # prod = [ 3.56*10^(-7), 1.0/18.0]
    # subS = [4,1]
    # prodS = [1,2]
    # ΔG0 = -47860.0
    ΔG = GFree(sub,subS,prod,prodS,T,ΔG0)
    println(ΔG)
    return(nothing)
end

@time main()
