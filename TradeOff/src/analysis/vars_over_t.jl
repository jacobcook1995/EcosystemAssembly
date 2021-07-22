# Script to find how variables change over time, which then saves them
using TradeOff
using JLD

function v_over_t()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 2
        error("insufficent inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    rps = 0
    ims = 0
    # Check that all arguments can be converted to integers
    try
        rps = parse(Int64,ARGS[1])
        ims = parse(Int64,ARGS[2])
    catch e
            error("need to provide 2 integers")
    end
    println("Compiled")
    flush(stdout)
    # Load in hardcoded simulation parameters
    Np, Nt, M, d = sim_paras()
    # Read in parameter file
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    # Load parameters
    ps = load(pfile,"ps")
    # List of pools already loaded in
    pls = []
    # Array of array to store pools
    pools = Array{Array{Microbe,1},1}(undef,1)
    # Counter for number of reactions
    NoR = 0
    # Loop over number of repeats
    for i = 1:rps
        # Load in relevant output file
        ofile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/Run$(i)Data$(ims)Ims.jld"
        if ~isfile(ofile)
            error("$(ims) immigrations run $(rN) is missing an output file")
        end
        # Load in microbe data, and immigration times
        T = load(ofile,"T")
        traj = load(ofile,"traj")
        micd = load(ofile,"micd")
        its = load(ofile,"its")
        # Use to construct full trajectory C
        C = merge_data(ps,traj,T,micd,its)
        # Preallocate vector of microbes
        ms = Array{Microbe,1}(undef,length(micd))
        # Loop over and find each one
        for j = 1:length(micd)
            # check for case where pool hasn't already been loaded in
            if micd[j].PID ∉ pls
                # Add new pool ID in
                pls = cat(pls,micd[j].PID,dims=1)
                # Find name of pool
                file = "Pools/ID=$(micd[j].PID)N=$(Nt)M=$(ps.M)d=$(d).jld"
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
        # Preallocate containers to store number of survivors with time
        svt = Array{Int64,1}(undef,length(T))
        tsvt = Array{Int64,1}(undef,length(T))
        sbs = Array{Int64,1}(undef,length(T))
        Rs = Array{Int64,2}(undef,NoR,length(T))
        via_R = Array{Int64,2}(undef,NoR,length(T))
        ηs = zeros(length(T))
        via_η = zeros(length(T))
        ωs = zeros(length(T))
        via_ω = zeros(length(T))
        fr_ΔG = zeros(length(T))
        ηs_R = zeros(NoR,length(T))
        ωs_R = zeros(NoR,length(T))
        # Save total number of strains
        numS = length(micd)
        # Loop over all time points
        for j = 1:length(T)
            # Find indices of surviving strains
            inds = findall(x->x>1e-5,C[j,1:numS])
            # Save number of surviving strains at each time point
            svt[j] = length(inds)
            # Find indices of "viable" strains
            vinds = findall(x->x>1e5,C[j,1:numS])
            # Save number of "viable" strains
            tsvt[j] = length(vinds)
            # Then also number of substrates
            sbs[j] = count(x->x>1e-12,C[j,(numS+1):(numS+ps.M)])
            # Loop over number of reactions
            for k = 1:NoR
                # Count number of strains with reaction for each case
                Rs[k,j] = count(x->x==k,ms[inds].↦:R)
                via_R[k,j] = count(x->x==k,ms[vinds].↦:R)
            end
            # Find (weighted) total eta value, and ω value
            for k = 1:length(inds)
                ηs[j] += sum(ms[inds[k]].η.*ms[inds[k]].ϕP)
                ωs[j] += ms[inds[k]].ω
            end
            # Average over number of strains
            if svt[j] > 0
                ηs[j] /= svt[j]
                ωs[j] /= svt[j]
            end
            # Find (weighted) total eta value for viable strains
            for k = 1:length(vinds)
                via_η[j] += sum(ms[vinds[k]].η.*ms[vinds[k]].ϕP)
                via_ω[j] += ms[vinds[k]].ω
            end
            # Average over number of strains
            if tsvt[j] > 0
                via_η[j] /= tsvt[j]
                via_ω[j] /= tsvt[j]
            end
            # Break down eta and omega value by R
            for k = 1:length(vinds)
                # Find relevant reaction number
                l = ms[vinds[k]].R
                # Add contribution to relevant total
                ηs_R[l,j] += sum(ms[vinds[k]].η.*ms[vinds[k]].ϕP)
                ωs_R[l,j] += ms[vinds[k]].ω
            end
            # Now weight by number of strains with each type of reaction
            for k = 1:NoR
                if via_R[k,j] > 0
                    ηs_R[k,j] /= via_R[k,j]
                    ωs_R[k,j] /= via_R[k,j]
                end
            end
            # Find fraction of free energy transduced
            for k = 1:length(vinds)
                for l = 1:ms[vinds[k]].R
                    # Find relevant reaction
                    r = ps.reacs[ms[vinds[k]].Reacs[l]]
                    # Then calculate frac
                    fr_ΔG[j] += ms[vinds[k]].η[l].*ms[vinds[k]].ϕP[l]*ΔGATP/(-r.ΔG0)
                end
            end
            # Average over number of viable strains
            if tsvt[j] > 0
                fr_ΔG[j] /= tsvt[j]
            end
        end
        # Now just save the relevant data
        jldopen("Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/AvRun$(i)Data$(ims)Ims.jld","w") do file
            # Save full timecourse
            write(file,"T",T)
            # Save reaction data
            write(file,"Rs",Rs)
            write(file,"via_R",via_R)
            write(file,"ηs_R",ηs_R)
            write(file,"ωs_R",ωs_R)
            # Save the other quantities
            write(file,"svt",svt)
            write(file,"tsvt",tsvt)
            write(file,"sbs",sbs)
            write(file,"ηs",ηs)
            write(file,"via_η",via_η)
            write(file,"ωs",ωs)
            write(file,"via_ω",via_ω)
            write(file,"fr_ΔG",fr_ΔG)
            # Finally save final time to help with benchmarking
            write(file,"Tf",T[end])
        end
        println("Run $i analysed")
        flush(stdout)
    end
    return(nothing)
end


@time v_over_t()
