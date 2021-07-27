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
    Np, Nt, M, d = sim_paras()
    # Number of steps to calculate stats for
    NumS = 2500
    # Read in parameter file
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    # Load parameters
    ps = load(pfile,"ps")
    # Container to store final times
    Tfs = zeros(rps)
    # Counter for number of reactions
    NoR = 0
    # Loop over number of repeats
    for i = 1:rps
        # Load in relevant output file
        vfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/AvRun$(i)Data$(ims)Ims.jld"
        if ~isfile(vfile)
            error("$(ims) immigrations run $(rN) is missing a variables file")
        end
        # Just want to save final times for now
        Tfs[i] = load(vfile,"Tf")
        # Save number of reactions from the first run
        if i == 1
            Rs = load(vfile,"Rs")
            NoR = size(Rs,1)
        end
    end
    # Use maximum final time to set time value
    times = collect(range(0.0,maximum(Tfs),length=NumS))
    # Preallocate relevant containers
    no_sims = zeros(length(times))
    cmb_svt = zeros(rps,length(times))
    cmb_tsvt = zeros(rps,length(times))
    cmb_sbs = zeros(rps,length(times))
    cmb_Rs = zeros(rps,NoR,length(times))
    cmb_via_R = zeros(rps,NoR,length(times))
    cmb_ηs_R = zeros(rps,NoR,length(times))
    cmb_ωs_R = zeros(rps,NoR,length(times))
    cmb_ηs = zeros(rps,length(times))
    cmb_via_η = zeros(rps,length(times))
    cmb_ωs = zeros(rps,length(times))
    cmb_via_ω = zeros(rps,length(times))
    cmb_fr_ΔG = zeros(rps,length(times))
    # Loop over number of trajectories (to minimise the number of reads in)
    for i = 1:rps
        # Load in relevant output file
        vfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/AvRun$(i)Data$(ims)Ims.jld"
        if ~isfile(vfile)
            error("$(ims) immigrations run $(rN) is missing a variables file")
        end
        T = load(vfile,"T")
        svt = load(vfile,"svt")
        tsvt = load(vfile,"tsvt")
        sbs = load(vfile,"sbs")
        Rs = load(vfile,"Rs")
        via_R = load(vfile,"via_R")
        ηs_R = load(vfile,"ηs_R")
        ωs_R = load(vfile,"ωs_R")
        ηs = load(vfile,"ηs")
        via_η = load(vfile,"via_η")
        ωs = load(vfile,"ωs")
        via_ω = load(vfile,"via_ω")
        fr_ΔG = load(vfile,"fr_ΔG")
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
                cmb_ηs[i,cnt] = ηs[Tind]*(T1x)/Tg + ηs[Tind-1]*(Tx2)/Tg
                cmb_via_η[i,cnt] = via_η[Tind]*(T1x)/Tg + via_η[Tind-1]*(Tx2)/Tg
                cmb_ωs[i,cnt] = ωs[Tind]*(T1x)/Tg + ωs[Tind-1]*(Tx2)/Tg
                cmb_via_ω[i,cnt] = via_ω[Tind]*(T1x)/Tg + via_ω[Tind-1]*(Tx2)/Tg
                cmb_fr_ΔG[i,cnt] = fr_ΔG[Tind]*(T1x)/Tg + fr_ΔG[Tind-1]*(Tx2)/Tg
                for j = 1:NoR
                    cmb_Rs[i,j,cnt] = Rs[j,Tind]*(T1x)/Tg + Rs[j,Tind-1]*(Tx2)/Tg
                    cmb_via_R[i,j,cnt] = via_R[j,Tind]*(T1x)/Tg + via_R[j,Tind-1]*(Tx2)/Tg
                    cmb_ηs_R[i,j,cnt] = ηs_R[j,Tind]*(T1x)/Tg + ηs_R[j,Tind-1]*(Tx2)/Tg
                    cmb_ωs_R[i,j,cnt] = ωs_R[j,Tind]*(T1x)/Tg + ωs_R[j,Tind-1]*(Tx2)/Tg
                end
            else
                # In the one case just add the value at time = 0
                cmb_svt[i,cnt] = svt[Tind]
                cmb_tsvt[i,cnt] = tsvt[Tind]
                cmb_sbs[i,cnt] = sbs[Tind]
                cmb_ηs[i,cnt] = ηs[Tind]
                cmb_via_η[i,cnt] = via_η[Tind]
                cmb_ωs[i,cnt] = ωs[Tind]
                cmb_via_ω[i,cnt] = via_ω[Tind]
                cmb_fr_ΔG[i,cnt] = fr_ΔG[Tind]
                for j = 1:NoR
                    cmb_Rs[i,j,cnt] = Rs[j,Tind]
                    cmb_via_R[i,j,cnt] = via_R[j,Tind]
                    cmb_ηs_R[i,j,cnt] = ηs_R[j,Tind]
                    cmb_ωs_R[i,j,cnt] = ωs_R[j,Tind]
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
    tot_via_R = dropdims(sum(cmb_via_R,dims=1),dims=1)
    tot_ηs_R = dropdims(sum(cmb_ηs_R,dims=1),dims=1)
    tot_ωs_R = dropdims(sum(cmb_ωs_R,dims=1),dims=1)
    tot_ηs = dropdims(sum(cmb_ηs,dims=1),dims=1)
    tot_via_η = dropdims(sum(cmb_via_η,dims=1),dims=1)
    tot_ωs = dropdims(sum(cmb_ωs,dims=1),dims=1)
    tot_via_ω = dropdims(sum(cmb_via_ω,dims=1),dims=1)
    tot_fr_ΔG = dropdims(sum(cmb_fr_ΔG,dims=1),dims=1)
    # Now calculate means
    mn_svt = tot_svt./no_sims
    mn_tsvt = tot_tsvt./no_sims
    mn_sbs = tot_sbs./no_sims
    mn_ηs = tot_ηs./no_sims
    mn_ωs = tot_ωs./no_sims
    # Preallocate viable strains data containers
    mn_via_ω = zeros(length(times))
    mn_via_η = zeros(length(times))
    mn_fr_ΔG = zeros(length(times))
    mn_via_R = zeros(NoR,length(times))
    # Loop over times to find number of viable strains
    for i = 1:length(times)
        # count number of cases with no viable strains
        nvs = count(x->x==0.0,cmb_tsvt[:,i])
        # Check if additional non viable strains have been found
        if nvs > rps - no_sims[i]
            vi_sim = rps - nvs
        else
            vi_sim = no_sims[i]
        end
        # Use to calculate accurate means
        mn_via_ω[i] = tot_via_ω[i]/(vi_sim)
        mn_via_η[i] = tot_via_η[i]/(vi_sim)
        mn_fr_ΔG[i] = tot_fr_ΔG[i]/(vi_sim)
        mn_via_R[:,i] = tot_via_R[:,i]./(vi_sim)
    end
    # 2D array has to be preallocated
    mn_Rs = zeros(NoR,length(times))
    mn_ηs_R = zeros(NoR,length(times))
    mn_ωs_R = zeros(NoR,length(times))
    for i = 1:NoR
        mn_Rs[i,:] = tot_Rs[i,:]./no_sims
        for j = 1:length(times)
            # count number of instances where reaction isn't present
            nvR = count(x->x==0.0,cmb_via_R[:,i,j])
            # Check if additional non viable strains have been found
            if nvR > rps - no_sims[i]
                vi_R = rps - nvR
            else
                vi_R = no_sims[i]
            end
            # Use number of strains with reaction to calculate the mean
            mn_ηs_R[i,j] = tot_ηs_R[i,j]/vi_R
            mn_ωs_R[i,j] = tot_ωs_R[i,j]/vi_R
        end
    end
    println("Means found")
    # Preallocate containers for the standard deviations
    sd_svt = zeros(size(mn_svt))
    sd_tsvt = zeros(size(mn_tsvt))
    sd_sbs = zeros(size(mn_sbs))
    sd_Rs = zeros(size(mn_Rs))
    sd_via_R = zeros(size(mn_via_R))
    sd_ηs_R = zeros(size(mn_ηs_R))
    sd_ωs_R = zeros(size(mn_ωs_R))
    sd_ηs = zeros(size(mn_ηs))
    sd_via_η = zeros(size(mn_via_η))
    sd_ωs = zeros(size(mn_ωs))
    sd_via_ω = zeros(size(mn_via_ω))
    sd_fr_ΔG = zeros(size(mn_fr_ΔG))
    # Loop over times
    for i = 1:length(times)
        # Find indices of still progressing trajectories
        inds = (Tfs .>= times[i])
        # Find indices of still progressing trajectories with one or more viable strains
        vinds = (Tfs .>= times[i]) .& (cmb_tsvt[:,i] .> 0.0)
        # calculate value to divide viable cases by
        no_via = sum(vinds)
        # Calculate standard deviations
        sd_svt[i] = sqrt(sum((cmb_svt[inds,i] .- mn_svt[i]).^2)/(no_sims[i] - 1))
        sd_tsvt[i] = sqrt(sum((cmb_tsvt[inds,i] .- mn_tsvt[i]).^2)/(no_sims[i] - 1))
        sd_sbs[i] = sqrt(sum((cmb_sbs[inds,i] .- mn_sbs[i]).^2)/(no_sims[i] - 1))
        sd_ηs[i] = sqrt(sum((cmb_ηs[inds,i] .- mn_ηs[i]).^2)/(no_sims[i] - 1))
        sd_ωs[i] = sqrt(sum((cmb_ωs[inds,i] .- mn_ωs[i]).^2)/(no_sims[i] - 1))
        # These should be calculated just for viable strains
        if no_via > 1
            sd_via_η[i] = sqrt(sum((cmb_via_η[vinds,i] .- mn_via_η[i]).^2)/(no_via - 1))
            sd_via_ω[i] = sqrt(sum((cmb_via_ω[vinds,i] .- mn_via_ω[i]).^2)/(no_via - 1))
            sd_fr_ΔG[i] = sqrt(sum((cmb_fr_ΔG[vinds,i] .- mn_fr_ΔG[i]).^2)/(no_via - 1))
            for j = 1:NoR
                sd_via_R[j,i] = sqrt(sum((cmb_via_R[vinds,j,i] .- mn_via_R[j,i]).^2)/(no_via - 1))
            end
        else
            sd_via_η[i] = NaN
            sd_via_ω[i] = NaN
            sd_fr_ΔG[i] = NaN
            sd_via_R[:,i] .= NaN
        end
        # Calculate standard deviations for reactions
        for j = 1:NoR
            sd_Rs[j,i] = sqrt(sum((cmb_Rs[inds,j,i] .- mn_Rs[j,i]).^2)/(no_sims[i] - 1))
            # Find indices of where reactions exist
            rinds = (Tfs .>= times[i]) .& (cmb_via_R[:,j,i] .> 0.0)
            # calculate value to divide viable cases by
            no_rs = sum(rinds)
            # Use only these in the reaction calculation
            if no_rs > 1
                sd_ηs_R[j,i] = sqrt(sum((cmb_ηs_R[rinds,j,i] .- mn_ηs_R[j,i]).^2)/(no_rs - 1))
                sd_ωs_R[j,i] = sqrt(sum((cmb_ωs_R[rinds,j,i] .- mn_ωs_R[j,i]).^2)/(no_rs - 1))
            else
                sd_ηs_R[j,i] = NaN
                sd_ωs_R[j,i] = NaN
            end
        end
    end
    # Now want to save means and standard deviations
    jldopen("Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/RunStats$(ims)Ims.jld","w") do file
        # Save times
        write(file,"times",times)
        # Save number of continuing trajectories
        write(file,"no_sims",no_sims)
        # Save averages
        write(file,"mn_svt",mn_svt)
        write(file,"mn_tsvt",mn_tsvt)
        write(file,"mn_sbs",mn_sbs)
        write(file,"mn_Rs",mn_Rs)
        write(file,"mn_via_R",mn_via_R)
        write(file,"mn_ηs_R",mn_ηs_R)
        write(file,"mn_ωs_R",mn_ωs_R)
        write(file,"mn_ηs",mn_ηs)
        write(file,"mn_via_η",mn_via_η)
        write(file,"mn_ωs",mn_ωs)
        write(file,"mn_via_ω",mn_via_ω)
        write(file,"mn_fr_ΔG",mn_fr_ΔG)
        # Save standard deviations
        write(file,"sd_svt",sd_svt)
        write(file,"sd_tsvt",sd_tsvt)
        write(file,"sd_sbs",sd_sbs)
        write(file,"sd_Rs",sd_Rs)
        write(file,"sd_via_R",sd_via_R)
        write(file,"sd_ηs_R",sd_ηs_R)
        write(file,"sd_ωs_R",sd_ωs_R)
        write(file,"sd_ηs",sd_ηs)
        write(file,"sd_via_η",sd_via_η)
        write(file,"sd_ωs",sd_ωs)
        write(file,"sd_via_ω",sd_via_ω)
        write(file,"sd_fr_ΔG",sd_fr_ΔG)
    end
    return(nothing)
end

@time trjstats()
