using Syntrophy
# This is a script to write my testing code into
# Anything reusable should be moved into Syntrophy.jl as a seperate callable function

function main()
    # Need to find a suitable testing example tomorrow
    T = 298.0 # Standard temp
    sub = [ 2.00*10^(-4)]
    prod = [ 3.00*10^(-6) ]
    subS = [1]
    prodS = [1]
    ΔG0 = 7549.6429
    ΔG = GFree(sub,subS,prod,prodS,T,ΔG0)
    println(ΔG/4184)
    return(nothing)
end

@time main()
