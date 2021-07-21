using TradeOff

# Function to generate a species pool
function make_spool()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 4
        error("Insufficent inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    Nt = 0
    M = 0
    Rl = 0
    Ru = 0
    # Check that all arguments can be converted to integers
    try
        Nt = parse(Int64,ARGS[1])
        M = parse(Int64,ARGS[2])
        Rl = parse(Int64,ARGS[3])
        Ru = parse(Int64,ARGS[4])
    catch e
            error("need to provide 4 integers")
    end
    # Check that simulation type is valid
    if Rl < 1
        error("lower bound on the number of reactions must be greater than 1")
    end
    # Check that simulation type is valid
    if Ru < Rl
        error("upper bound on the number of reactions can't be smaller than the lower")
    end
    # Check that number of strains is greater than 0
    if Nt < 1
        error("number of strains should be greater than zero")
    end
    # Check number of metabolites is greater than 1
    if M < 2
        error("model requires at least two metabolites")
    end
    # Make desired vector of reactions
    Rs = collect(Rl:Ru)
    # HARDCODING THESE IN FOR NOW, WORK MORE CAREFULLY ON THESE LATER
    d = 6e-5
    μrange = 5e6*(M/25)
    mratio = 1e-2 # Product to substrate ratio for equilbrium
    # Finally generate the new pool
    new_pool(Nt,M,Rs,d,μrange,mratio)
    return(nothing)
end

@time make_spool()
