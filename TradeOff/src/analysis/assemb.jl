# A script to assemble and save communities
using TradeOff
using JLD
using Glob

# Function to assemble specfic communities
function assemble()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 1
        error("Insufficent inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    rps = 0
    # Check that all arguments can be converted to integers
    try
        rps = parse(Int64,ARGS[1])
    catch e
            error("need to provide 1 integer")
    end
    # Starting run assumed to be 1
    Rs = 1
    # Check if run to start from has been provided
    if length(ARGS) > 1
        # Check this argument is an integer
        try
            Rs = parse(Int64,ARGS[2])
        catch e
            error("intial run number must be integer")
        end
    end
    # Check that number of strains is greater than 0
    if rps < 1
        error("need to do at least 1 simulation")
    end
    if Rs > rps
        error("starting run can't be higher than final run")
    end
    # Now read in hard coded simulation parameters
    Np, Rls, Rus, Nt, M = sim_paras()
    # Now start actual script
    println("Compiled and input read in!")
    flush(stdout)
    # Use formula to find how many reactions number of metabolites implies
    O = 2*M - 3
    # Preallocate container for filenames
    pls = fill("",Np)
    # Loop over number of required pools
    for i = 1:Np
        # Find all pools satisfying the condition
        flnms = glob("Pools/ID=*N=$(Nt)M=$(M)Reacs$(Rls[i])-$(Rus[i]).jld")
        # Loop over valid filenames
        for j = 1:length(flnms)
            # Save first that hasn't already been used
            if flnms[j] ∉ pls
                pls[i] = flnms[j]
            end
        end
    end
    # Save the reaction set for the first file as a point of comparison
    rs = load(pls[1],"reacs")
    # Check that all pools match this
    for i = 2:Np
        rst = load(pls[i],"reacs")
        if rst ≠ rs
            error("pool $i uses different reaction set")
        end
    end
    # Initial ribosome fraction is taken from ATP fits I did a while ago
    ϕR0 = 0.128
    # Fairly arbitary inital conditions
    pop = 1000.0
    conc = 0.0
    as = 1e5
    ϕs = ϕR0
    # Starting with 10 strains for now
    Ni = 10
    # Mean immigration time assumed to be 1*10^5 seconds
    mT = 1e5
    # Small number of immigration events for inital testing
    ims = 250
    # Rate of additional immigrants
    λIm = 0.5
    # Make parameter set
    ps = initialise(M,O)
    # Check that reaction set is identical to sets the pool was generated with
    if ps.reacs ≠ rs
        error("simulation reaction set does not match pool reaction set")
    end
    # Check if directory exists and if not make it
    if ~isdir("Output/$(Np)Pools$(M)Metabolites$(Nt)Species")
        mkdir("Output/$(Np)Pools$(M)Metabolites$(Nt)Species")
    end
    # Save this parameter set
    jldopen("Output/$(Np)Pools$(M)Metabolites$(Nt)Species/Paras$(ims)Ims.jld","w") do file
        write(file,"ps",ps)
    end
    # ONLY LOADING ONE POOL AT THE MOMENT, THIS PROBABLY HAS TO CHANGE
    mpl = load(pls[1],"mics")
    # Now loop over the number of repeats
    for i = Rs:rps
        # Print that the new run has been started
        println("Run $i started!")
        flush(stdout)
        # Find starting time
        ti = time()
        # Then run the simulation
        C, T, micd, its = full_simulate(ps,pop,conc,as,ϕs,mpl,Ni,mT,ims,λIm)
        # And then print time elapsed
        tf = time()
        println("Time elapsed on run $i: $(tf-ti) s")
        # Now just save the relevant data
        jldopen("Output/$(Np)Pools$(M)Metabolites$(Nt)Species/Run$(i)Data$(ims)Ims.jld","w") do file
            # Save full set of microbe data
            write(file,"micd",micd)
            # Save extinction times
            write(file,"its",its)
            # Save time data and dynamics data
            write(file,"T",T)
            write(file,"C",C)
        end
        # Print to show that run has been successfully completed
        println("Run $i completed and saved!")
        flush(stdout)
    end
    return(nothing)
end

# Hard code simulation parameters into this function
function sim_paras()
    # Set the hardcoded variables here
    Np = 1
    Rls = [1]
    Rus = [5]
    Nt = 1000
    M = 25
    return(Np,Rls,Rus,Nt,M)
end

@time assemble()
