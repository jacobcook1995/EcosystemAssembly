using Syntrophy
# This is a script to write my testing code into
# Anything reusable should be moved into Syntrophy.jl as a seperate callable function

function main()
    greet()
    return(nothing)
end

@time main()
