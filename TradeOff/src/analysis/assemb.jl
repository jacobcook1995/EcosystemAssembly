# A script to assemble and save communities
using TradeOff
using JLD
using Glob

# WILL PROBABLY NEED TO CHANGE SIGNIFICANTLY WHEN I CHANGE THE ASSEMBLY PROCEDURE

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
    # NEED TO COMPLETELY RETHINK THE SOMETHING DYNAMICS HERE
    # FUNDAMENTALLY NEED TO FIGURE OUT HOW CALLBACKS WORK BEFORE I CAN PROCEED HERE
    # Set time long enough for dynamics to equilbrate
    Tmax = 1e8
    # THESE PROBABLY CAN GO BUT I NEED TO FIGURE OUT WHAT TO REPLACE THEM WITH
    # Initial ribosome fraction is taken from my ATP fits
    ϕR0 = 0.128
    # Fairly arbitary inital conditions
    pop = 1.0
    conc = 0.0
    as = 1e5
    ϕs = ϕR0
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
    jldopen("Output/$(Np)Pools$(M)Metabolites$(Nt)Species/Paras.jld","w") do file
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
        C, T = full_simulate(ps,Tmax,pop,conc,as,ϕs,mpl)
        # And then print time elapsed
        tf = time()
        println("Time elapsed on run $i: $(tf-ti) s")
        flush(stdout)
        return(nothing)
        # THIS STUFF IS PROBABLY REDUNDANT NOW, BUT NEED SOMETHING TO REPLACE IT WITH
        # Establish which microbes are extinct
        ext = (C[end,1:N] .== 0.0)
        # Preallocate vector to store extinct microbes
        ded = Array{Microbe,1}(undef,sum(ext))
        # Loop over and store microbes in the vector
        k = 0
        for j = 1:length(ext)
            if ext[j] == 1
                k += 1
                ded[k] = ps.mics[j]
            end
        end
        # Remove extinct strains from parameter set
        ps = extinction(ps,ext)
        # Preallocate final concentrations (etc) for output
        out = Array{Float64,1}(undef,3*ps.N+M)
        # Store final metabolite concentrations
        out[ps.N+1:ps.N+M] = C[end,N+1:N+M]
        # Now sub in data for not extinct microbes
        k = 0
        for j = 1:length(ext)
            if ext[j] != 1
                k += 1
                # Population
                out[k] = C[end,j]
                # Energy
                out[M+ps.N+k] = C[end,M+N+j]
                # Fraction
                out[M+2*ps.N+k] = C[end,M+2*N+j]
            end
        end
        # Save extinct strains
        jldopen("Output/ExtinctReacs$(Rl)-$(Ru)Run$(i)Ns$(N).jld","w") do file
            write(file,"ded",ded)
        end
        # and the full output
        jldopen("Output/OutputReacs$(Rl)-$(Ru)Run$(i)Ns$(N).jld","w") do file
            # Save final output
            write(file,"out",out)
            # # Save time data and dynamics data
            write(file,"T",T)
            write(file,"C",C[1:end,1:end])
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
