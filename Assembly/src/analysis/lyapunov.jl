# A script to run Lyapunov analysis of our inhibition model.
using Assembly
using JLD

# function to set up and run a lyapunov analysis of a simple case
function lya()
    # only a small number of strains and metabolites for this analysis
    N = 2
    M = 2
    O = 1
    mR = 1.0
    sdR = 0.0
    mq = 1.0
    sdq = 0.1
    mK  = 0.1
    sdK = 0.01
    mk = 10.0
    sdk = 1.0
    # Now make the parameter set, can just use the normal function for now
    ps = initialise(N,M,O,mR,sdR,mq,sdq,mK,sdK,mk,sdk)
    # save this parameter set
    jldopen("Temp/Paras/psL.jld","w") do file
        write(file,"ps",ps)
    end
    return(nothing)
end

@time lya()
