# Script to find and analyse the interaction strength and type of strains in my simulation
using Assembly
using JLD
using SymPy

# function to quantify the interactions
function quantify_ints()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 4
        error("Insufficent inputs provided (looking for 4)")
    end
    # Preallocate the variables I want to extract from the input
    Rl = 0
    Ru = 0
    syn = true
    rps = 0
    # Check that all arguments can be converted to integers
    try
        Rl = parse(Int64,ARGS[1])
        Ru = parse(Int64,ARGS[2])
        syn = parse(Bool,ARGS[3])
        rps = parse(Int64,ARGS[4])
    catch e
           error("need to provide 3 integers and a bool")
    end
    # Check that simulation type is valid
    if Rl < 1
        error("lower bound has to be at least one reaction")
    end
    if Ru < Rl
        error("upper bound can't be lower than the lower bound")
    end
    # Check that number of simulations is greater than 0
    if rps < 1
        error("number of repeats can't be less than 1")
    end
    println("Compiled!")
    # For testing purposes just want to look at a specific simulation
    # For 1-5 true case simulation 4 has just one survivor
    # In this case simulation 6 has two suvivors
    for i = 6 #1:rps
        # Read in relevant files
        pfile = "Data/$(Rl)-$(Ru)$(syn)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        ofile = "Data/$(Rl)-$(Ru)$(syn)/RedOutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        if ~isfile(ofile)
            error("run $(i) is missing an output file")
        end
        efile = "Data/$(Rl)-$(Ru)$(syn)/RedExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        if ~isfile(efile)
            error("run $(i) is missing an extinct file")
        end
        # Basically just loading everything out as I'm not sure what I'll need
        ps = load(pfile,"ps")
        C = load(ofile,"C")
        T = load(ofile,"T")
        out = load(ofile,"out")
        inf_out = load(ofile,"inf_out")
        ded = load(efile,"ded")
        # Preallocate vector of forces
        F = Array{Sym,1}(undef,3*ps.N+ps.M)
        # Now find vector of forces
        F = Force(ps,F)
    end
    return(nothing)
end


@time quantify_ints()
