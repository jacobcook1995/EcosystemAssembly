# File to calculate stats for the variables across trajectories
using TradeOff
using JLD

# Function to read in variables with time and calculate stats over time
function trjstats()
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
    # Token to insert into filenames
    tk = ""
    # Overwritten for no immigration case
    if sim_type == 5
        tk = "NoImm"
    end
    # Load in hardcoded simulation parameters
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Number of steps to calculate stats for
    NumS = 2500
    # Read in parameter file
    pfile = "Output/$(tk)$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations run is missing a parameter file")
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
        vfile = "Output/$(tk)$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/AvRun$(i)Data$(ims)Ims.jld"
        if ~isfile(vfile)
            error("$(ims) immigrations run $(i) is missing a variables file")
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
    no_via = zeros(length(times))
    no_rs = zeros(NoR,length(times))
    cmb_svt = zeros(rps,length(times))
    cmb_tsvt = zeros(rps,length(times))
    cmb_pop = zeros(rps,length(times))
    cmb_via_bm = zeros(rps,length(times))
    cmb_shD = zeros(rps,length(times))
    cmb_sbs = zeros(rps,length(times))
    cmb_Rs = zeros(rps,NoR,length(times))
    cmb_via_R = zeros(rps,NoR,length(times))
    cmb_ηs_R = zeros(rps,NoR,length(times))
    cmb_ωs_R = zeros(rps,NoR,length(times))
    cmb_kc_R = zeros(rps,NoR,length(times))
    cmb_KS_R = zeros(rps,NoR,length(times))
    cmb_kr_R = zeros(rps,NoR,length(times))
    cmb_ηs = zeros(rps,length(times))
    cmb_via_η = zeros(rps,length(times))
    cmb_via_η_bw = zeros(rps,length(times))
    cmb_ωs = zeros(rps,length(times))
    cmb_via_ω = zeros(rps,length(times))
    cmb_via_ω_bw = zeros(rps,length(times))
    cmb_fr_ΔG = zeros(rps,length(times))
    cmb_fr_ΔG_bw = zeros(rps,length(times))
    cmb_via_a = zeros(rps,length(times))
    cmb_via_ϕR = zeros(rps,length(times))
    cmb_kcs = zeros(rps,length(times))
    cmb_KSs = zeros(rps,length(times))
    cmb_krs = zeros(rps,length(times))
    cmb_av_steps = zeros(rps,length(times))
    cmb_av_steps_bw = zeros(rps,length(times))
    cmb_η_stp = zeros(rps,length(times),M-1)
    cmb_fr_ΔG_stp = zeros(rps,length(times),M-1)
    cmb_ϕP_stp = zeros(rps,length(times),M-1)
    all_fin_ϕRs = Float64[]
    all_l_sb = zeros(rps)
    # Loop over number of trajectories (to minimise the number of reads in)
    for i = 1:rps
        # Load in relevant output file
        vfile = "Output/$(tk)$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/AvRun$(i)Data$(ims)Ims.jld"
        if ~isfile(vfile)
            error("$(ims) immigrations run $(rN) is missing a variables file")
        end
        T = load(vfile,"T")
        svt = load(vfile,"svt")
        tsvt = load(vfile,"tsvt")
        pop = load(vfile,"pop")
        via_bm = load(vfile,"via_bm")
        shD = load(vfile,"shD")
        sbs = load(vfile,"sbs")
        Rs = load(vfile,"Rs")
        via_R = load(vfile,"via_R")
        ηs_R = load(vfile,"ηs_R")
        ωs_R = load(vfile,"ωs_R")
        kc_R = load(vfile,"kc_R")
        KS_R = load(vfile,"KS_R")
        kr_R = load(vfile,"kr_R")
        ηs = load(vfile,"ηs")
        via_η = load(vfile,"via_η")
        via_η_bw = load(vfile,"via_η_bw")
        ωs = load(vfile,"ωs")
        via_ω = load(vfile,"via_ω")
        via_ω_bw = load(vfile,"via_ω_bw")
        fr_ΔG = load(vfile,"fr_ΔG")
        fr_ΔG_bw = load(vfile,"fr_ΔG_bw")
        via_a = load(vfile,"via_a")
        via_ϕR = load(vfile,"via_ϕR")
        kcs = load(vfile,"kcs")
        KSs = load(vfile,"KSs")
        krs = load(vfile,"krs")
        av_steps = load(vfile,"av_steps")
        av_steps_bw = load(vfile,"av_steps_bw")
        η_stp = load(vfile,"η_stp")
        fr_ΔG_stp = load(vfile,"fr_ΔG_stp")
        ϕP_stp = load(vfile,"ϕP_stp")
        fin_ϕR = load(vfile,"fin_ϕR")
        all_fin_ϕRs = cat(all_fin_ϕRs,fin_ϕR,dims=1)
        all_l_sb[i] = load(vfile,"l_sb")
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
            # Only increment the viable counter if there are viable strains at time point
            if tsvt[Tind] != 0
                no_via[cnt] += 1
            end
            # Loop over reactions to find number of viable reactions
            for j = 1:NoR
                # Check if counter should be incremented
                if via_R[j,Tind] .> 0.0
                    no_rs[j,cnt] += 1
                end
            end
            # Skip averaging if previous point is missing
            if Tind > 1 && tsvt[Tind-1] != 0
                # Calculate relevant time gaps
                Tg = (T[Tind]-T[Tind-1])
                T1x = times[cnt]-T[Tind-1]
                Tx2 = T[Tind]-times[cnt]
                # And use to find appropriate averages
                cmb_svt[i,cnt] = svt[Tind]*(T1x)/Tg + svt[Tind-1]*(Tx2)/Tg
                cmb_tsvt[i,cnt] = tsvt[Tind]*(T1x)/Tg + tsvt[Tind-1]*(Tx2)/Tg
                cmb_pop[i,cnt] = pop[Tind]*(T1x)/Tg + pop[Tind-1]*(Tx2)/Tg
                cmb_via_bm[i,cnt] = via_bm[Tind]*(T1x)/Tg + via_bm[Tind-1]*(Tx2)/Tg
                cmb_shD[i,cnt] = shD[Tind]*(T1x)/Tg + shD[Tind-1]*(Tx2)/Tg
                cmb_sbs[i,cnt] = sbs[Tind]*(T1x)/Tg + sbs[Tind-1]*(Tx2)/Tg
                cmb_ηs[i,cnt] = ηs[Tind]*(T1x)/Tg + ηs[Tind-1]*(Tx2)/Tg
                cmb_via_η[i,cnt] = via_η[Tind]*(T1x)/Tg + via_η[Tind-1]*(Tx2)/Tg
                cmb_via_η_bw[i,cnt] = via_η_bw[Tind]*(T1x)/Tg + via_η_bw[Tind-1]*(Tx2)/Tg
                cmb_ωs[i,cnt] = ωs[Tind]*(T1x)/Tg + ωs[Tind-1]*(Tx2)/Tg
                cmb_via_ω[i,cnt] = via_ω[Tind]*(T1x)/Tg + via_ω[Tind-1]*(Tx2)/Tg
                cmb_via_ω_bw[i,cnt] = via_ω_bw[Tind]*(T1x)/Tg + via_ω_bw[Tind-1]*(Tx2)/Tg
                cmb_fr_ΔG[i,cnt] = fr_ΔG[Tind]*(T1x)/Tg + fr_ΔG[Tind-1]*(Tx2)/Tg
                cmb_fr_ΔG_bw[i,cnt] = fr_ΔG_bw[Tind]*(T1x)/Tg + fr_ΔG_bw[Tind-1]*(Tx2)/Tg
                cmb_via_a[i,cnt] = via_a[Tind]*(T1x)/Tg + via_a[Tind-1]*(Tx2)/Tg
                cmb_via_ϕR[i,cnt] = via_ϕR[Tind]*(T1x)/Tg + via_ϕR[Tind-1]*(Tx2)/Tg
                cmb_kcs[i,cnt] = kcs[Tind]*(T1x)/Tg + kcs[Tind-1]*(Tx2)/Tg
                cmb_KSs[i,cnt] = KSs[Tind]*(T1x)/Tg + KSs[Tind-1]*(Tx2)/Tg
                cmb_krs[i,cnt] = krs[Tind]*(T1x)/Tg + krs[Tind-1]*(Tx2)/Tg
                cmb_av_steps[i,cnt] = av_steps[Tind]*(T1x)/Tg + av_steps[Tind-1]*(Tx2)/Tg
                cmb_av_steps_bw[i,cnt] = av_steps_bw[Tind]*(T1x)/Tg + av_steps_bw[Tind-1]*(Tx2)/Tg
                cmb_η_stp[i,cnt,:] = η_stp[Tind,:]*(T1x)/Tg .+ η_stp[Tind-1,:]*(Tx2)/Tg
                cmb_fr_ΔG_stp[i,cnt,:] = fr_ΔG_stp[Tind,:]*(T1x)/Tg .+ fr_ΔG_stp[Tind-1,:]*(Tx2)/Tg
                cmb_ϕP_stp[i,cnt,:] = ϕP_stp[Tind,:]*(T1x)/Tg .+ ϕP_stp[Tind-1,:]*(Tx2)/Tg
                for j = 1:NoR
                    cmb_Rs[i,j,cnt] = Rs[j,Tind]*(T1x)/Tg + Rs[j,Tind-1]*(Tx2)/Tg
                    cmb_via_R[i,j,cnt] = via_R[j,Tind]*(T1x)/Tg + via_R[j,Tind-1]*(Tx2)/Tg
                    cmb_ηs_R[i,j,cnt] = ηs_R[j,Tind]*(T1x)/Tg + ηs_R[j,Tind-1]*(Tx2)/Tg
                    cmb_ωs_R[i,j,cnt] = ωs_R[j,Tind]*(T1x)/Tg + ωs_R[j,Tind-1]*(Tx2)/Tg
                    cmb_kc_R[i,j,cnt] = kc_R[j,Tind]*(T1x)/Tg + kc_R[j,Tind-1]*(Tx2)/Tg
                    cmb_KS_R[i,j,cnt] = KS_R[j,Tind]*(T1x)/Tg + KS_R[j,Tind-1]*(Tx2)/Tg
                    cmb_kr_R[i,j,cnt] = kr_R[j,Tind]*(T1x)/Tg + kr_R[j,Tind-1]*(Tx2)/Tg
                end
            else
                # In the one case just add the value at time = 0
                cmb_svt[i,cnt] = svt[Tind]
                cmb_tsvt[i,cnt] = tsvt[Tind]
                cmb_pop[i,cnt] = pop[Tind]
                cmb_via_bm[i,cnt] = via_bm[Tind]
                cmb_shD[i,cnt] = shD[Tind]
                cmb_sbs[i,cnt] = sbs[Tind]
                cmb_ηs[i,cnt] = ηs[Tind]
                cmb_via_η[i,cnt] = via_η[Tind]
                cmb_via_η_bw[i,cnt] = via_η_bw[Tind]
                cmb_ωs[i,cnt] = ωs[Tind]
                cmb_via_ω[i,cnt] = via_ω[Tind]
                cmb_via_ω_bw[i,cnt] = via_ω_bw[Tind]
                cmb_fr_ΔG[i,cnt] = fr_ΔG[Tind]
                cmb_fr_ΔG_bw[i,cnt] = fr_ΔG_bw[Tind]
                cmb_via_a[i,cnt] = via_a[Tind]
                cmb_via_ϕR[i,cnt] = via_ϕR[Tind]
                cmb_kcs[i,cnt] = kcs[Tind]
                cmb_KSs[i,cnt] = KSs[Tind]
                cmb_krs[i,cnt] = krs[Tind]
                cmb_av_steps[i,cnt] = av_steps[Tind]
                cmb_av_steps_bw[i,cnt] = av_steps_bw[Tind]
                cmb_η_stp[i,cnt,:] = η_stp[Tind,:]
                cmb_fr_ΔG_stp[i,cnt,:] = fr_ΔG_stp[Tind,:]
                cmb_ϕP_stp[i,cnt,:] = ϕP_stp[Tind,:]
                for j = 1:NoR
                    cmb_Rs[i,j,cnt] = Rs[j,Tind]
                    cmb_via_R[i,j,cnt] = via_R[j,Tind]
                    cmb_ηs_R[i,j,cnt] = ηs_R[j,Tind]
                    cmb_ωs_R[i,j,cnt] = ωs_R[j,Tind]
                    cmb_kc_R[i,j,cnt] = kc_R[j,Tind]
                    cmb_KS_R[i,j,cnt] = KS_R[j,Tind]
                    cmb_kr_R[i,j,cnt] = kr_R[j,Tind]
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
    tot_pop = dropdims(sum(cmb_pop,dims=1),dims=1)
    tot_via_bm = dropdims(sum(cmb_via_bm,dims=1),dims=1)
    tot_shD = dropdims(sum(cmb_shD,dims=1),dims=1)
    tot_sbs = dropdims(sum(cmb_sbs,dims=1),dims=1)
    tot_Rs = dropdims(sum(cmb_Rs,dims=1),dims=1)
    tot_via_R = dropdims(sum(cmb_via_R,dims=1),dims=1)
    tot_ηs_R = dropdims(sum(cmb_ηs_R,dims=1),dims=1)
    tot_ωs_R = dropdims(sum(cmb_ωs_R,dims=1),dims=1)
    tot_kc_R = dropdims(sum(cmb_kc_R,dims=1),dims=1)
    tot_KS_R = dropdims(sum(cmb_KS_R,dims=1),dims=1)
    tot_kr_R = dropdims(sum(cmb_kr_R,dims=1),dims=1)
    tot_ηs = dropdims(sum(cmb_ηs,dims=1),dims=1)
    tot_via_η = dropdims(sum(cmb_via_η,dims=1),dims=1)
    tot_via_η_bw = dropdims(sum(cmb_via_η_bw,dims=1),dims=1)
    tot_ωs = dropdims(sum(cmb_ωs,dims=1),dims=1)
    tot_via_ω = dropdims(sum(cmb_via_ω,dims=1),dims=1)
    tot_via_ω_bw = dropdims(sum(cmb_via_ω_bw,dims=1),dims=1)
    tot_fr_ΔG = dropdims(sum(cmb_fr_ΔG,dims=1),dims=1)
    tot_fr_ΔG_bw = dropdims(sum(cmb_fr_ΔG_bw,dims=1),dims=1)
    tot_via_a = dropdims(sum(cmb_via_a,dims=1),dims=1)
    tot_via_ϕR = dropdims(sum(cmb_via_ϕR,dims=1),dims=1)
    tot_kcs = dropdims(sum(cmb_kcs,dims=1),dims=1)
    tot_KSs = dropdims(sum(cmb_KSs,dims=1),dims=1)
    tot_krs = dropdims(sum(cmb_krs,dims=1),dims=1)
    tot_av_steps = dropdims(sum(cmb_av_steps,dims=1),dims=1)
    tot_av_steps_bw = dropdims(sum(cmb_av_steps_bw,dims=1),dims=1)
    tot_η_stp = dropdims(sum(cmb_η_stp,dims=1),dims=1)
    tot_fr_ΔG_stp = dropdims(sum(cmb_fr_ΔG_stp,dims=1),dims=1)
    tot_ϕP_stp = dropdims(sum(cmb_ϕP_stp,dims=1),dims=1)
    # Now calculate means
    mn_svt = tot_svt./no_sims
    mn_tsvt = tot_tsvt./no_sims
    mn_pop = tot_pop./no_sims
    mn_shD = tot_shD./no_sims
    mn_sbs = tot_sbs./no_sims
    mn_ηs = tot_ηs./no_sims
    mn_ωs = tot_ωs./no_sims
    # Calculate means for viable case
    mn_via_bm = tot_via_bm./(no_via)
    mn_via_ω = tot_via_ω./(no_via)
    mn_via_ω_bw = tot_via_ω_bw./(no_via)
    mn_via_η = tot_via_η./(no_via)
    mn_via_η_bw = tot_via_η_bw./(no_via)
    mn_via_a = tot_via_a./(no_via)
    mn_via_ϕR = tot_via_ϕR./(no_via)
    mn_fr_ΔG = tot_fr_ΔG./(no_via)
    mn_fr_ΔG_bw = tot_fr_ΔG_bw./(no_via)
    mn_kcs = tot_kcs./(no_via)
    mn_KSs = tot_KSs./(no_via)
    mn_krs = tot_krs./(no_via)
    mn_av_steps = tot_av_steps./(no_via)
    mn_av_steps_bw = tot_av_steps_bw./(no_via)
    mn_η_stp = tot_η_stp./(no_via)
    mn_fr_ΔG_stp = tot_fr_ΔG_stp./(no_via)
    mn_ϕP_stp = tot_ϕP_stp./(no_via)
    # Preallocate 2D array
    mn_via_R = zeros(NoR,length(times))
    for i = 1:NoR
        mn_via_R[i,:] = tot_via_R[i,:]./(no_via)
    end
    # 2D arrays have to be preallocated
    mn_Rs = zeros(NoR,length(times))
    mn_ηs_R = zeros(NoR,length(times))
    mn_ωs_R = zeros(NoR,length(times))
    mn_kc_R = zeros(NoR,length(times))
    mn_KS_R = zeros(NoR,length(times))
    mn_kr_R = zeros(NoR,length(times))
    for i = 1:NoR
        mn_Rs[i,:] = tot_Rs[i,:]./no_sims
        # Use number of strains with reaction to calculate the mean
        mn_ηs_R[i,:] = tot_ηs_R[i,:]./no_rs[i,:]
        mn_ωs_R[i,:] = tot_ωs_R[i,:]./no_rs[i,:]
        mn_kc_R[i,:] = tot_kc_R[i,:]./no_rs[i,:]
        mn_KS_R[i,:] = tot_KS_R[i,:]./no_rs[i,:]
        mn_kr_R[i,:] = tot_kr_R[i,:]./no_rs[i,:]
    end
    println("Means found")
    # Preallocate containers for the standard deviations
    sd_svt = zeros(size(mn_svt))
    sd_tsvt = zeros(size(mn_tsvt))
    sd_pop = zeros(size(mn_pop))
    sd_via_bm = zeros(size(mn_via_bm))
    sd_shD = zeros(size(mn_shD))
    sd_sbs = zeros(size(mn_sbs))
    sd_Rs = zeros(size(mn_Rs))
    sd_via_R = zeros(size(mn_via_R))
    sd_ηs_R = zeros(size(mn_ηs_R))
    sd_ωs_R = zeros(size(mn_ωs_R))
    sd_kc_R = zeros(size(mn_kc_R))
    sd_KS_R = zeros(size(mn_KS_R))
    sd_kr_R = zeros(size(mn_kr_R))
    sd_ηs = zeros(size(mn_ηs))
    sd_via_η = zeros(size(mn_via_η))
    sd_via_η_bw = zeros(size(mn_via_η_bw))
    sd_ωs = zeros(size(mn_ωs))
    sd_via_ω = zeros(size(mn_via_ω))
    sd_via_ω_bw = zeros(size(mn_via_ω_bw))
    sd_fr_ΔG = zeros(size(mn_fr_ΔG))
    sd_fr_ΔG_bw = zeros(size(mn_fr_ΔG_bw))
    sd_via_a = zeros(size(mn_via_a))
    sd_via_ϕR = zeros(size(mn_via_ϕR))
    sd_kcs = zeros(size(mn_kcs))
    sd_KSs = zeros(size(mn_KSs))
    sd_krs = zeros(size(mn_krs))
    sd_av_steps = zeros(size(mn_av_steps))
    sd_av_steps_bw = zeros(size(mn_av_steps_bw))
    sd_η_stp = zeros(size(mn_η_stp))
    sd_fr_ΔG_stp = zeros(size(mn_fr_ΔG_stp))
    sd_ϕP_stp = zeros(size(mn_ϕP_stp))
    # Loop over times
    for i = 1:length(times)
        # Find indices of still progressing trajectories
        inds = (Tfs .>= times[i])
        # Find indices of still progressing trajectories with one or more viable strains
        vinds = (Tfs .>= times[i]) .& (cmb_tsvt[:,i] .> 0.0)
        # Calculate standard deviations
        sd_svt[i] = sqrt(sum((cmb_svt[inds,i] .- mn_svt[i]).^2)/(no_sims[i] - 1))
        sd_tsvt[i] = sqrt(sum((cmb_tsvt[inds,i] .- mn_tsvt[i]).^2)/(no_sims[i] - 1))
        sd_pop[i] = sqrt(sum((cmb_pop[inds,i] .- mn_pop[i]).^2)/(no_sims[i] - 1))
        sd_shD[i] = sqrt(sum((cmb_shD[inds,i] .- mn_shD[i]).^2)/(no_sims[i] - 1))
        sd_sbs[i] = sqrt(sum((cmb_sbs[inds,i] .- mn_sbs[i]).^2)/(no_sims[i] - 1))
        sd_ηs[i] = sqrt(sum((cmb_ηs[inds,i] .- mn_ηs[i]).^2)/(no_sims[i] - 1))
        sd_ωs[i] = sqrt(sum((cmb_ωs[inds,i] .- mn_ωs[i]).^2)/(no_sims[i] - 1))
        # These should be calculated just for viable strains
        if no_via[i] > 1
            sd_via_bm[i] = sqrt(sum((cmb_via_bm[inds,i] .- mn_via_bm[i]).^2)/(no_via[i] - 1))
            sd_via_η[i] = sqrt(sum((cmb_via_η[vinds,i] .- mn_via_η[i]).^2)/(no_via[i] - 1))
            sd_via_η_bw[i] = sqrt(sum((cmb_via_η_bw[vinds,i] .- mn_via_η_bw[i]).^2)/(no_via[i] - 1))
            sd_via_ω[i] = sqrt(sum((cmb_via_ω[vinds,i] .- mn_via_ω[i]).^2)/(no_via[i] - 1))
            sd_via_ω_bw[i] = sqrt(sum((cmb_via_ω_bw[vinds,i] .- mn_via_ω_bw[i]).^2)/(no_via[i] - 1))
            sd_fr_ΔG[i] = sqrt(sum((cmb_fr_ΔG[vinds,i] .- mn_fr_ΔG[i]).^2)/(no_via[i] - 1))
            sd_fr_ΔG_bw[i] = sqrt(sum((cmb_fr_ΔG_bw[vinds,i] .- mn_fr_ΔG_bw[i]).^2)/(no_via[i] - 1))
            sd_via_a[i] = sqrt(sum((cmb_via_a[vinds,i] .- mn_via_a[i]).^2)/(no_via[i] - 1))
            sd_via_ϕR[i] = sqrt(sum((cmb_via_ϕR[vinds,i] .- mn_via_ϕR[i]).^2)/(no_via[i] - 1))
            sd_kcs[i] = sqrt(sum((cmb_kcs[vinds,i] .- mn_kcs[i]).^2)/(no_via[i] - 1))
            sd_KSs[i] = sqrt(sum((cmb_KSs[vinds,i] .- mn_KSs[i]).^2)/(no_via[i] - 1))
            sd_krs[i] = sqrt(sum((cmb_krs[vinds,i] .- mn_krs[i]).^2)/(no_via[i] - 1))
            sd_av_steps[i] = sqrt(sum((cmb_av_steps[vinds,i] .- mn_av_steps[i]).^2)/(no_via[i] - 1))
            sd_av_steps_bw[i] = sqrt(sum((cmb_av_steps_bw[vinds,i] .- mn_av_steps_bw[i]).^2)/(no_via[i] - 1))
            for j = 1:(M-1)
                sd_η_stp[i,j] = sqrt(sum((cmb_η_stp[vinds,i,j] .- mn_η_stp[i,j]).^2)/(no_via[i] - 1))
                sd_fr_ΔG_stp[i,j] = sqrt(sum((cmb_fr_ΔG_stp[vinds,i,j] .- mn_fr_ΔG_stp[i,j]).^2)/(no_via[i] - 1))
                sd_ϕP_stp[i,j] = sqrt(sum((cmb_ϕP_stp[vinds,i,j] .- mn_ϕP_stp[i,j]).^2)/(no_via[i] - 1))
            end
            for j = 1:NoR
                sd_via_R[j,i] = sqrt(sum((cmb_via_R[vinds,j,i] .- mn_via_R[j,i]).^2)/(no_via[i] - 1))
            end
        else
            sd_via_bm[i] = NaN
            sd_via_η[i] = NaN
            sd_via_η_bw[i] = NaN
            sd_via_ω[i] = NaN
            sd_via_ω_bw[i] = NaN
            sd_fr_ΔG[i] = NaN
            sd_fr_ΔG_bw[i] = NaN
            sd_via_a[i] = NaN
            sd_via_ϕR[i] = NaN
            sd_kcs[i] = NaN
            sd_KSs[i] = NaN
            sd_krs[i] = NaN
            sd_av_steps[i] = NaN
            sd_av_steps_bw[i] = NaN
            sd_η_stp[i,:] .= NaN
            sd_fr_ΔG_stp[i,:] .= NaN
            sd_ϕP_stp[i,:] .= NaN
            sd_via_R[:,i] .= NaN
        end
        # Calculate standard deviations for reactions
        for j = 1:NoR
            sd_Rs[j,i] = sqrt(sum((cmb_Rs[inds,j,i] .- mn_Rs[j,i]).^2)/(no_sims[i] - 1))
            # Find indices of where reactions exist
            rinds = (Tfs .>= times[i]) .& (cmb_via_R[:,j,i] .> 0.0)
            # Use only these in the reaction calculation
            if no_rs[j,i] > 1
                sd_ηs_R[j,i] = sqrt(sum((cmb_ηs_R[rinds,j,i] .- mn_ηs_R[j,i]).^2)/(no_rs[j,i] - 1))
                sd_ωs_R[j,i] = sqrt(sum((cmb_ωs_R[rinds,j,i] .- mn_ωs_R[j,i]).^2)/(no_rs[j,i] - 1))
                sd_kc_R[j,i] = sqrt(sum((cmb_kc_R[rinds,j,i] .- mn_kc_R[j,i]).^2)/(no_rs[j,i] - 1))
                sd_KS_R[j,i] = sqrt(sum((cmb_KS_R[rinds,j,i] .- mn_KS_R[j,i]).^2)/(no_rs[j,i] - 1))
                sd_kr_R[j,i] = sqrt(sum((cmb_kr_R[rinds,j,i] .- mn_kr_R[j,i]).^2)/(no_rs[j,i] - 1))
            else
                sd_ηs_R[j,i] = NaN
                sd_ωs_R[j,i] = NaN
                sd_kc_R[j,i] = NaN
                sd_KS_R[j,i] = NaN
                sd_kr_R[j,i] = NaN
            end
        end
    end
    # Now want to save means and standard deviations
    jldopen("Output/$(tk)$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/RunStats$(ims)Ims.jld","w") do file
        # Save times
        write(file,"times",times)
        # Save number of continuing trajectories
        write(file,"no_sims",no_sims)
        write(file,"no_via",no_via)
        write(file,"no_rs",no_rs)
        # Save averages
        write(file,"mn_svt",mn_svt)
        write(file,"mn_tsvt",mn_tsvt)
        write(file,"mn_pop",mn_pop)
        write(file,"mn_via_bm",mn_via_bm)
        write(file,"mn_shD",mn_shD)
        write(file,"mn_sbs",mn_sbs)
        write(file,"mn_Rs",mn_Rs)
        write(file,"mn_via_R",mn_via_R)
        write(file,"mn_ηs_R",mn_ηs_R)
        write(file,"mn_ωs_R",mn_ωs_R)
        write(file,"mn_kc_R",mn_kc_R)
        write(file,"mn_KS_R",mn_KS_R)
        write(file,"mn_kr_R",mn_kr_R)
        write(file,"mn_ηs",mn_ηs)
        write(file,"mn_via_η",mn_via_η)
        write(file,"mn_via_η_bw",mn_via_η_bw)
        write(file,"mn_ωs",mn_ωs)
        write(file,"mn_via_ω",mn_via_ω)
        write(file,"mn_via_ω_bw",mn_via_ω_bw)
        write(file,"mn_fr_ΔG",mn_fr_ΔG)
        write(file,"mn_fr_ΔG_bw",mn_fr_ΔG_bw)
        write(file,"mn_via_a",mn_via_a)
        write(file,"mn_via_ϕR",mn_via_ϕR)
        write(file,"mn_kcs",mn_kcs)
        write(file,"mn_KSs",mn_KSs)
        write(file,"mn_krs",mn_krs)
        write(file,"mn_av_steps",mn_av_steps)
        write(file,"mn_av_steps_bw",mn_av_steps_bw)
        write(file,"mn_η_stp",mn_η_stp)
        write(file,"mn_fr_ΔG_stp",mn_fr_ΔG_stp)
        write(file,"mn_ϕP_stp",mn_ϕP_stp)
        # Save standard deviations
        write(file,"sd_svt",sd_svt)
        write(file,"sd_tsvt",sd_tsvt)
        write(file,"sd_pop",sd_pop)
        write(file,"sd_via_bm",sd_via_bm)
        write(file,"sd_shD",sd_shD)
        write(file,"sd_sbs",sd_sbs)
        write(file,"sd_Rs",sd_Rs)
        write(file,"sd_via_R",sd_via_R)
        write(file,"sd_ηs_R",sd_ηs_R)
        write(file,"sd_ωs_R",sd_ωs_R)
        write(file,"sd_kc_R",sd_kc_R)
        write(file,"sd_KS_R",sd_KS_R)
        write(file,"sd_kr_R",sd_kr_R)
        write(file,"sd_ηs",sd_ηs)
        write(file,"sd_via_η",sd_via_η)
        write(file,"sd_via_η_bw",sd_via_η_bw)
        write(file,"sd_ωs",sd_ωs)
        write(file,"sd_via_ω",sd_via_ω)
        write(file,"sd_via_ω_bw",sd_via_ω_bw)
        write(file,"sd_fr_ΔG",sd_fr_ΔG)
        write(file,"sd_fr_ΔG_bw",sd_fr_ΔG_bw)
        write(file,"sd_via_a",sd_via_a)
        write(file,"sd_via_ϕR",sd_via_ϕR)
        write(file,"sd_kcs",sd_kcs)
        write(file,"sd_KSs",sd_KSs)
        write(file,"sd_krs",sd_krs)
        write(file,"sd_av_steps",sd_av_steps)
        write(file,"sd_av_steps_bw",sd_av_steps_bw)
        write(file,"sd_η_stp",sd_η_stp)
        write(file,"sd_fr_ΔG_stp",sd_fr_ΔG_stp)
        write(file,"sd_ϕP_stp",sd_ϕP_stp)
        # Write all of the ϕ values out
        write(file,"all_fin_ϕRs",all_fin_ϕRs)
        write(file,"all_l_sb",all_l_sb)
    end
    return(nothing)
end

# Function to read in snapshot data and calculate stats over time
function snpstats()
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
        error("$(ims) immigrations run is missing a parameter file")
    end
    # Load parameters
    ps = load(pfile,"ps")
    # Load in 1st output file
    sfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/SnapData$(ims)Ims.jld"
    if ~isfile(sfile)
        error("$(ims) immigrations run 1 is missing a snapshot data file")
    end
    # Load in snapshot times
    times = load(sfile,"times")
    ns = load(sfile,"ns")
    gs = load(sfile,"gs")
    stb = load(sfile,"stb")
    inc = load(sfile,"inc")
    dec = load(sfile,"dec")
    st_r = load(sfile,"st_r")
    # Connstruct totals here
    tot_stb = dropdims(sum(stb,dims=2),dims=2)
    tot_inc = dropdims(sum(inc,dims=2),dims=2)
    tot_dec = dropdims(sum(dec,dims=2),dims=2)
    # Now calculate means
    mn_stb = tot_stb./st_r
    mn_inc = tot_inc./st_r
    mn_dec = tot_dec./st_r
    println("Means found")
    # Preallocate containers for the standard deviations
    sd_stb = zeros(size(mn_stb))
    sd_inc = zeros(size(mn_inc))
    sd_dec = zeros(size(mn_dec))
    # Loop over times
    for i = 1:length(times)-1
        # Find indices of still progressing trajectories
        inds = (ns[i,:] .!== 0.0)
        # Calculate standard deviations
        sd_stb[i] = sqrt(sum((stb[i,inds] .- mn_stb[i]).^2)/(st_r[i] - 1))
        sd_inc[i] = sqrt(sum((inc[i,inds] .- mn_inc[i]).^2)/(st_r[i] - 1))
        sd_dec[i] = sqrt(sum((dec[i,inds] .- mn_dec[i]).^2)/(st_r[i] - 1))
    end
    # Find total numbers of new strains and total number that grows
    tot_ns = dropdims(sum(ns,dims=2),dims=2)
    tot_gs = dropdims(sum(gs,dims=2),dims=2)
    # Use to find growth probability
    gp = zeros(size(tot_gs))
    for i = 1:length(gp)
        # Check that value isn't zero
        if tot_ns[i] != 0
            gp[i] = tot_gs[i]/tot_ns[i]
        else
            gp[i] = 0.0
        end
    end
    # Now just save the relevant data
    jldopen("Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/SnapDataStats$(ims)Ims.jld","w") do file
        # Save times of snapshots
        write(file,"times",times)
        # Save growth probabilities
        write(file,"gp",gp)
        # Save means
        write(file,"mn_stb",mn_stb)
        write(file,"mn_inc",mn_inc)
        write(file,"mn_dec",mn_dec)
        # Save standard deviations
        write(file,"sd_stb",sd_stb)
        write(file,"sd_inc",sd_inc)
        write(file,"sd_dec",sd_dec)
        # Save number of simulations
        write(file,"st_r",st_r)
    end
    return(nothing)
end

@time trjstats()
