# Script to calculate features of the trajectories at specific snapshot times
using TradeOff
using JLD

# Function to calculate relevant values at partciular time snap shots
function snp_shot()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 3
        error("insufficent inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    rps = 0
    ims = 0
    sim_type = 0
    # Check that all arguments can be converted to integers
    try
        rps = parse(Int64,ARGS[1])
        ims = parse(Int64,ARGS[2])
        sim_type = parse(Int64,ARGS[3])
    catch e
            error("need to provide 3 integers")
    end
    println("Compiled")
    flush(stdout)
    # Load in hardcoded simulation parameters
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Read in parameter file
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations is missing a parameter file")
    end
    # Load parameters
    ps = load(pfile,"ps")
    # Preallocate vector to store final times
    Tfs = zeros(rps)
    # Loop over number of repeats to find maximum times
    for i = 1:rps
        # Load in relevant output file
        ofile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Run$(i)Data$(ims)Ims.jld"
        if ~isfile(ofile)
            error("$(ims) immigrations run $(i) is missing an output file")
        end
        # Load in microbe data, and immigration times
        T = load(ofile,"T")
        # Store final T value
        Tfs[i] = T[end]
    end
    # Then find maximum final time
    Tmax = maximum(Tfs)
    # Number of steps to calculate stats for
    NumS = 2500
    # Define snap shot times based on this maximum time
    snps = collect(range(0.0,Tmax,length=NumS))
    # Also save the time step for later use
    t_step = snps[2] - snps[1]
    # Preallocate data to save
    ns = zeros(length(snps)-1,rps)
    gs = zeros(length(snps)-1,rps)
    # List of pools already loaded in
    pls = []
    # Array of array to store pools
    pools = Array{Array{Microbe,1},1}(undef,1)
    # Counter for number of reactions
    NoR = 0
    # Loop over number of repeats
    for i = 1:rps
        # Load in relevant output file
        ofile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Run$(i)Data$(ims)Ims.jld"
        if ~isfile(ofile)
            error("$(ims) immigrations run $(i) is missing an output file")
        end
        # Load in microbe data, and immigration times
        T = load(ofile,"T")
        traj = load(ofile,"traj")
        micd = load(ofile,"micd")
        its = load(ofile,"its")
        # Use to construct full trajectory C
        C = merge_data(ps,traj,T,micd,its)
        # Find and save initial population value for this run
        Ni = C[1,1]
        # Preallocate vector of microbes
        ms = Array{Microbe,1}(undef,length(micd))
        # Loop over and find each one
        for j = 1:length(micd)
            # check for case where pool hasn't already been loaded in
            if micd[j].PID ∉ pls
                # Add new pool ID in
                pls = cat(pls,micd[j].PID,dims=1)
                # Find name of pool
                file = "Pools/ID=$(micd[j].PID)N=$(Nt)M=$(ps.M)d=$(d)u=$(μrange).jld"
                # Check if this is the first pool
                if length(pls) == 1
                    # If so save the pool
                    pools[1] = load(file,"mics")
                    # Find number of reactions based on this
                    NoR = maximum(pools[1].↦:R)
                else
                    # Otherwise just cat it on existing vector
                    pool = load(file,"mics")
                    pools = cat(pools,pool,dims=1)
                    # Find maximum number of reactions for this pool
                    NoRt = maximum(pools[1].↦:R)
                    # Save if higher than old number of reactions
                    NoR = max(NoR,NoRt)
                end
            end
            # Find correct pool to read from
            ind = findfirst(x->x==micd[j].PID,pls)
            # Use this index to find and save the correct microbe
            ms[j] = (pools[ind])[micd[j].MID]
        end
        # loop over time snapshots
        for j = 1:length(snps)-1
            # Find indices of the time point before and after the snapshot point
            ind1 = findfirst(x->x>=snps[j],T)
            ind2 = findfirst(x->x>=snps[j+1],T)
            # Check for new species that entered the system in this time window
            migs = findall(x->snps[j]<=x<snps[j+1],micd.↦:ImT)
            # Count number of new immigrants
            ns[j,i] = length(migs)
            # Also need to calculate the number that grow (over the snapshot period)
            for k = 1:length(migs)
                # Check that population has increased from the initial value
                if C[ind2,migs[k]] > Ni
                    gs[j,i] += 1
                end
            end
        end
        # Now just save the relevant data
        jldopen("Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/SnapRun$(i)Data$(ims)Ims.jld","w") do file
            # Save times of snapshots
            write(file,"times",snps)
            # Save whatever I generate here
            write(file,"ns",ns)
            write(file,"gs",gs)
            # Finally save final time to help with benchmarking
            write(file,"Tf",T[end])
        end
        println("Run $i analysed")
        flush(stdout)
        return(nothing)
    end
    return(nothing)
end

@time snp_shot()
