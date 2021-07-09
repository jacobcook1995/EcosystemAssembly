# File to calculate stats for the variables across trajectories
using TradeOff
using JLD

# Function to read in variables with time and calculate stats over time
function trjstats()
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
    Np, Rls, Rus, Nt, M = sim_paras()
    # Save number of reactions
    NoR = Rus[1] - Rls[1] + 1
    # Number of steps to calculate stats for
    NumS = 2500
    # Read in parameter file
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Species/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    # Load parameters
    ps = load(pfile,"ps")
    # Container to store final times
    Tfs = zeros(rps)
    # Loop over number of repeats
    for i = 1:rps
        # Load in relevant output file
        vfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Species/AvRun$(i)Data$(ims)Ims.jld"
        if ~isfile(vfile)
            error("$(ims) immigrations run $(rN) is missing a variables file")
        end
        # Just want to save final times for now
        Tfs[i] = load(vfile,"Tf")
    end
    # Use maximum final time to set time value
    times = collect(range(0.0,maximum(Tfs),length=NumS))
    # Preallocate relevant containers
    no_sims = zeros(length(times))
    cmb_svt = zeros(rps,length(times))
    cmb_tsvt = zeros(rps,length(times))
    cmb_sbs = zeros(rps,length(times))
    cmb_Rs = zeros(rps,NoR,length(times))
    # Loop over number of trajectories (to minimise the number of reads in)
    for i = 1:rps
        # Load in relevant output file
        vfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Species/AvRun$(i)Data$(ims)Ims.jld"
        if ~isfile(vfile)
            error("$(ims) immigrations run $(rN) is missing a variables file")
        end
        T = load(vfile,"T")
        svt = load(vfile,"svt")
        tsvt = load(vfile,"tsvt")
        sbs = load(vfile,"sbs")
        Rs = load(vfile,"Rs")
        # Bool to indicate end of the run
        Rn_end = false
        cnt = 0
        # Loop until end of run reached
        while Rn_end == false
            # increment counter
            cnt += 1
            # Still a trajectory here so increment that count by one
            no_sims[cnt] += 1
            # Find index of first point greater than or equal to time
            Tind = findfirst(x->x>=times[cnt],T)
            if Tind > 1
                # Calculate relevant time gaps
                Tg = (T[Tind]-T[Tind-1])
                T1x = times[cnt]-T[Tind-1]
                Tx2 = T[Tind]-times[cnt]
                # And use to find appropriate averages
                cmb_svt[i,cnt] = svt[Tind]*(T1x)/Tg + svt[Tind-1]*(Tx2)/Tg
                cmb_tsvt[i,cnt] = tsvt[Tind]*(T1x)/Tg + tsvt[Tind-1]*(Tx2)/Tg
                cmb_sbs[i,cnt] = sbs[Tind]*(T1x)/Tg + sbs[Tind-1]*(Tx2)/Tg
                for j = 1:NoR
                    cmb_Rs[i,j,cnt] = Rs[j,Tind]*(T1x)/Tg + Rs[j,Tind-1]*(Tx2)/Tg
                end
            else
                # In the one case just add the value at time = 0
                cmb_svt[i,cnt] = svt[Tind]
                cmb_tsvt[i,cnt] = tsvt[Tind]
                cmb_sbs[i,cnt] = sbs[Tind]
                for j = 1:NoR
                    cmb_Rs[i,j,cnt] = Rs[j,Tind]
                end
            end
            # Finally check if next time point is higher than final time for this trajectory
            if cnt >= length(times) || times[cnt+1] > Tfs[i]
                Rn_end = true
            end
        end
        println("Analysed trajectory $(i)")
        flush(stdout)
    end
    # Sum to find totals
    tot_svt = dropdims(sum(cmb_svt,dims=1),dims=1)
    tot_tsvt = dropdims(sum(cmb_tsvt,dims=1),dims=1)
    tot_sbs = dropdims(sum(cmb_sbs,dims=1),dims=1)
    tot_Rs = dropdims(sum(cmb_Rs,dims=1),dims=1)
    # Now calculate means
    mn_svt = tot_svt./no_sims
    mn_tsvt = tot_tsvt./no_sims
    mn_sbs = tot_sbs./no_sims
    # 2D array has to be preallocated
    mn_Rs = zeros(NoR,length(times))
    for i = 1:NoR
        mn_Rs[i,:] = tot_Rs[i,:]./no_sims
    end
    println("Means found")
    # Preallocate containers for the standard deviations
    sd_svt = zeros(size(mn_svt))
    sd_tsvt = zeros(size(mn_tsvt))
    sd_sbs = zeros(size(mn_sbs))
    sd_Rs = zeros(size(mn_Rs))
    # Loop over times
    for i = 1:length(times)
        # Find indices of still progressing trajectories
        inds = (Tfs .>= times[i])
        # Calculate standard deviations
        sd_svt[i] = sqrt(sum((cmb_svt[inds,i] .- mn_svt[i]).^2)/(no_sims[i] - 1))
        sd_tsvt[i] = sqrt(sum((cmb_tsvt[inds,i] .- mn_tsvt[i]).^2)/(no_sims[i] - 1))
        sd_sbs[i] = sqrt(sum((cmb_sbs[inds,i] .- mn_sbs[i]).^2)/(no_sims[i] - 1))
        for j = 1:NoR
            sd_Rs[j,i] = sqrt(sum((cmb_Rs[inds,j,i] .- mn_Rs[j,i]).^2)/(no_sims[i] - 1))
        end
    end
    # Now want to save means and standard deviations
    jldopen("Output/$(Np)Pools$(M)Metabolites$(Nt)Species/RunStats$(ims)Ims.jld","w") do file
        # Save times
        write(file,"times",times)
        # Save number of continuing trajectories
        write(file,"no_sims",no_sims)
        # Save averages
        write(file,"mn_svt",mn_svt)
        write(file,"mn_tsvt",mn_tsvt)
        write(file,"mn_sbs",mn_sbs)
        write(file,"mn_Rs",mn_Rs)
        # Save standard deviations
        write(file,"sd_svt",sd_svt)
        write(file,"sd_tsvt",sd_tsvt)
        write(file,"sd_sbs",sd_sbs)
        write(file,"sd_Rs",sd_Rs)
    end
    return(nothing)
end

@time trjstats()
