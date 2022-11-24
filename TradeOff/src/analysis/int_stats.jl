# File to calculate stats for the variables across trajectories
using TradeOff
using JLD

# Function to interpolate over a time series
function interpolate_time(ts::Array[Float64,1], Tg::Float64, T1x::Float64, T2x::Float64)
    return (ts[Tind] * (T1x) / Tg + ts[Tind-1] * (T2x) / Tg)

# Function to read in variables with time and calculate stats over time
function intstats()
    # Check that sufficient arguments have been provided
    if length(ARGS) < 2
        error("insufficient inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    rps = 0
    ims = 0
    # Check that all arguments can be converted to integers
    try
        rps = parse(Int64, ARGS[1])
        ims = parse(Int64, ARGS[2])
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
    ps = load(pfile, "ps")
    # Container to store final times
    Tfs = zeros(rps)
    # Counter for number of reactions
    NoR = 0
    # Loop over number of repeats
    for i = 1:rps
        # Load in relevant output file
        vfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/IntsRun$(i)Data$(ims)Ims.jld"
        if ~isfile(vfile)
            error("$(ims) immigrations run $(rN) is missing a variables file")
        end
        # Just want to save final times for now
        Tfs[i] = load(vfile, "Tf")
    end
    # Use maximum final time to set time value
    times = collect(range(0.0, maximum(Tfs), length = NumS))
    # Preallocate relevant containers
    no_sims = zeros(length(times))
    no_via = zeros(length(times))
    cmb_svt = zeros(rps, length(times))
    cmb_tsvt = zeros(rps, length(times))
    cmb_no_comp = zeros(rps, length(times))
    cmb_no_facl = zeros(rps, length(times))
    cmb_no_selff = zeros(rps, length(times))
    cmb_no_selfc = zeros(rps, length(times))
    cmb_via_no_comp = zeros(rps, length(times))
    cmb_via_no_facl = zeros(rps, length(times))
    cmb_via_no_selff = zeros(rps, length(times))
    cmb_via_no_selfc = zeros(rps, length(times))
    cmb_st_comp = zeros(rps, length(times))
    cmb_st_facl = zeros(rps, length(times))
    cmb_st_selfc = zeros(rps, length(times))
    cmb_st_selff = zeros(rps, length(times))
    # Loop over number of trajectories (to minimise the number of reads in)
    for i = 1:rps
        # Load in relevant output file
        vfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/IntsRun$(i)Data$(ims)Ims.jld"
        if ~isfile(vfile)
            error("$(ims) immigrations run $(rN) is missing a variables file")
        end
        T = load(vfile, "T")
        svt = load(vfile, "svt")
        tsvt = load(vfile, "tsvt")
        no_comp = load(vfile, "no_comp")
        no_facl = load(vfile, "no_facl")
        no_selff = load(vfile, "no_selff")
        no_selfc = load(vfile, "no_selfc")
        via_no_comp = load(vfile, "via_no_comp")
        via_no_facl = load(vfile, "via_no_facl")
        via_no_selff = load(vfile, "via_no_selff")
        via_no_selfc = load(vfile, "via_no_selfc")
        st_comp = load(vfile, "st_comp")
        st_facl = load(vfile, "st_facl")
        st_selfc = load(vfile, "st_selfc")
        st_selff = load(vfile, "st_selff")
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
            Tind = findfirst(x -> x >= times[cnt], T)
            # Only increment the viable counter if there are viable strains at time point
            if tsvt[Tind] != 0
                no_via[cnt] += 1
            end
            # Update the all but the first step in the same way
            if Tind > 1
                # Calculate relevant time gaps
                Tg = (T[Tind] - T[Tind-1])
                T1x = times[cnt] - T[Tind-1]
                T2x = T[Tind] - times[cnt]
                # And use to find appropriate averages
                cmb_svt[i, cnt] = interpolate_time(svt, Tg, T1x, T2x)
                cmb_tsvt[i, cnt] = interpolate_time(tsvt, Tg, T1x, T2x)
                cmb_no_comp[i, cnt] = interpolate_time(no_comp, Tg, T1x, T2x)
                cmb_no_facl[i, cnt] = interpolate_time(no_facl, Tg, T1x, T2x)
                cmb_no_selff[i, cnt] = interpolate_time(no_selff, Tg, T1x, T2x)
                cmb_no_selfc[i, cnt] = interpolate_time(no_selfc, Tg, T1x, T2x)
                cmb_via_no_comp[i, cnt] = interpolate_time(via_no_comp, Tg, T1x, T2x)
                cmb_via_no_facl[i, cnt] = interpolate_time(via_no_facl, Tg, T1x, T2x)
                cmb_via_no_selff[i, cnt] = interpolate_time(via_no_selff, Tg, T1x, T2x)
                cmb_via_no_selfc[i, cnt] = interpolate_time(via_no_selfc, Tg, T1x, T2x)
                cmb_st_comp[i, cnt] = interpolate_time(via_st_comp, Tg, T1x, T2x)
                cmb_st_facl[i, cnt] = interpolate_time(via_st_facl, Tg, T1x, T2x)
                cmb_st_selfc[i, cnt] = interpolate_time(via_st_selfc, Tg, T1x, T2x)
                cmb_st_selff[i, cnt] = interpolate_time(via_st_selff, Tg, T1x, T2x)
            else
                # In the one case just add the value at time = 0
                cmb_svt[i, cnt] = svt[Tind]
                cmb_tsvt[i, cnt] = tsvt[Tind]
                cmb_no_comp[i, cnt] = no_comp[Tind]
                cmb_no_facl[i, cnt] = no_facl[Tind]
                cmb_no_selff[i, cnt] = no_selff[Tind]
                cmb_no_selfc[i, cnt] = no_selfc[Tind]
                cmb_via_no_comp[i, cnt] = via_no_comp[Tind]
                cmb_via_no_facl[i, cnt] = via_no_facl[Tind]
                cmb_via_no_selff[i, cnt] = via_no_selff[Tind]
                cmb_via_no_selfc[i, cnt] = via_no_selfc[Tind]
                cmb_st_comp[i, cnt] = st_comp[Tind]
                cmb_st_facl[i, cnt] = st_facl[Tind]
                cmb_st_selfc[i, cnt] = st_selfc[Tind]
                cmb_st_selff[i, cnt] = st_selff[Tind]
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
    tot_svt = dropdims(sum(cmb_svt, dims = 1), dims = 1)
    tot_tsvt = dropdims(sum(cmb_tsvt, dims = 1), dims = 1)
    tot_no_comp = dropdims(sum(cmb_no_comp, dims = 1), dims = 1)
    tot_no_facl = dropdims(sum(cmb_no_facl, dims = 1), dims = 1)
    tot_no_selff = dropdims(sum(cmb_no_selff, dims = 1), dims = 1)
    tot_no_selfc = dropdims(sum(cmb_no_selfc, dims = 1), dims = 1)
    tot_via_no_comp = dropdims(sum(cmb_via_no_comp, dims = 1), dims = 1)
    tot_via_no_facl = dropdims(sum(cmb_via_no_facl, dims = 1), dims = 1)
    tot_via_no_selff = dropdims(sum(cmb_via_no_selff, dims = 1), dims = 1)
    tot_via_no_selfc = dropdims(sum(cmb_via_no_selfc, dims = 1), dims = 1)
    tot_st_comp = dropdims(sum(cmb_st_comp, dims = 1), dims = 1)
    tot_st_facl = dropdims(sum(cmb_st_facl, dims = 1), dims = 1)
    tot_st_selfc = dropdims(sum(cmb_st_selfc, dims = 1), dims = 1)
    tot_st_selff = dropdims(sum(cmb_st_selff, dims = 1), dims = 1)
    # Now calculate means
    mn_svt = tot_svt ./ no_sims
    mn_tsvt = tot_tsvt ./ no_sims
    mn_no_comp = tot_no_comp ./ no_sims
    mn_no_facl = tot_no_facl ./ no_sims
    mn_no_selff = tot_no_selff ./ no_sims
    mn_no_selfc = tot_no_selfc ./ no_sims
    mn_st_comp = tot_st_comp ./ no_sims
    mn_st_facl = tot_st_facl ./ no_sims
    mn_st_selfc = tot_st_selfc ./ no_sims
    mn_st_selff = tot_st_selff ./ no_sims
    # Calculate for the viable case
    mn_via_no_comp = tot_via_no_comp ./ no_via
    mn_via_no_facl = tot_via_no_facl ./ no_via
    mn_via_no_selfc = tot_via_no_selfc ./ no_via
    mn_via_no_selff = tot_via_no_selff ./ no_via
    println("Means found")
    # Preallocate containers for the standard deviations
    sd_svt = zeros(size(mn_svt))
    sd_tsvt = zeros(size(mn_tsvt))
    sd_no_comp = zeros(size(tot_no_comp))
    sd_no_facl = zeros(size(tot_no_facl))
    sd_no_selff = zeros(size(tot_no_selff))
    sd_no_selfc = zeros(size(tot_no_selfc))
    sd_via_no_comp = zeros(size(tot_via_no_comp))
    sd_via_no_facl = zeros(size(tot_via_no_facl))
    sd_via_no_selff = zeros(size(tot_via_no_selff))
    sd_via_no_selfc = zeros(size(tot_via_no_selfc))
    sd_st_comp = zeros(size(tot_st_comp))
    sd_st_facl = zeros(size(tot_st_facl))
    sd_st_selfc = zeros(size(tot_st_selfc))
    sd_st_selff = zeros(size(tot_st_selff))
    # Loop over times
    for i in eachindex(times)
        # Find indices of still progressing trajectories
        inds = (Tfs .>= times[i])
        # Find indices of still progressing trajectories with one or more viable strains
        vinds = (Tfs .>= times[i]) .& (cmb_tsvt[:, i] .> 0.0)
        # Calculate standard deviations
        sd_svt[i] = sqrt(sum((cmb_svt[inds, i] .- mn_svt[i]) .^ 2) / (no_sims[i] - 1))
        sd_tsvt[i] = sqrt(sum((cmb_tsvt[inds, i] .- mn_tsvt[i]) .^ 2) / (no_sims[i] - 1))
        sd_no_comp[i] =
            sqrt(sum((cmb_no_comp[inds, i] .- mn_no_comp[i]) .^ 2) / (no_sims[i] - 1))
        sd_no_facl[i] =
            sqrt(sum((cmb_no_facl[inds, i] .- mn_no_facl[i]) .^ 2) / (no_sims[i] - 1))
        sd_no_selff[i] =
            sqrt(sum((cmb_no_selff[inds, i] .- mn_no_selff[i]) .^ 2) / (no_sims[i] - 1))
        sd_no_selfc[i] =
            sqrt(sum((cmb_no_selfc[inds, i] .- mn_no_selfc[i]) .^ 2) / (no_sims[i] - 1))
        sd_st_comp[i] =
            sqrt(sum((cmb_st_comp[inds, i] .- mn_st_comp[i]) .^ 2) / (no_sims[i] - 1))
        sd_st_facl[i] =
            sqrt(sum((cmb_st_facl[inds, i] .- mn_st_facl[i]) .^ 2) / (no_sims[i] - 1))
        sd_st_selff[i] =
            sqrt(sum((cmb_st_selff[inds, i] .- mn_st_selff[i]) .^ 2) / (no_sims[i] - 1))
        sd_st_selfc[i] =
            sqrt(sum((cmb_st_selfc[inds, i] .- mn_st_selfc[i]) .^ 2) / (no_sims[i] - 1))
        # These should be calculated just for viable strains
        if no_via[i] > 1
            sd_via_no_comp[i] = sqrt(
                sum((cmb_via_no_comp[vinds, i] .- mn_via_no_comp[i]) .^ 2) /
                (no_via[i] - 1),
            )
            sd_via_no_facl[i] = sqrt(
                sum((cmb_via_no_facl[vinds, i] .- mn_via_no_facl[i]) .^ 2) /
                (no_via[i] - 1),
            )
            sd_via_no_selff[i] = sqrt(
                sum((cmb_via_no_selff[vinds, i] .- mn_via_no_selff[i]) .^ 2) /
                (no_via[i] - 1),
            )
            sd_via_no_selfc[i] = sqrt(
                sum((cmb_via_no_selfc[vinds, i] .- mn_via_no_selfc[i]) .^ 2) /
                (no_via[i] - 1),
            )
        else
            sd_via_no_comp[i] = NaN
            sd_via_no_facl[i] = NaN
            sd_via_no_selff[i] = NaN
            sd_via_no_selfc[i] = NaN
        end
    end
    # Now want to save means and standard deviations
    jldopen(
        "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/IntStats$(ims)Ims.jld",
        "w",
    ) do file
        # Save times
        write(file, "times", times)
        # Save number of continuing trajectories
        write(file, "no_sims", no_sims)
        write(file, "no_via", no_via)
        # Save averages
        write(file, "mn_svt", mn_svt)
        write(file, "mn_tsvt", mn_tsvt)
        write(file, "mn_no_comp", mn_no_comp)
        write(file, "mn_no_facl", mn_no_facl)
        write(file, "mn_no_selff", mn_no_selff)
        write(file, "mn_no_selfc", mn_no_selfc)
        write(file, "mn_via_no_comp", mn_via_no_comp)
        write(file, "mn_via_no_facl", mn_via_no_facl)
        write(file, "mn_via_no_selff", mn_via_no_selff)
        write(file, "mn_via_no_selfc", mn_via_no_selfc)
        write(file, "mn_st_comp", mn_st_comp)
        write(file, "mn_st_facl", mn_st_facl)
        write(file, "mn_st_selfc", mn_st_selfc)
        write(file, "mn_st_selff", mn_st_selff)
        # Save standard deviations
        write(file, "sd_svt", sd_svt)
        write(file, "sd_tsvt", sd_tsvt)
        write(file, "sd_no_comp", sd_no_comp)
        write(file, "sd_no_facl", sd_no_facl)
        write(file, "sd_no_selff", sd_no_selff)
        write(file, "sd_no_selfc", sd_no_selfc)
        write(file, "sd_via_no_comp", sd_via_no_comp)
        write(file, "sd_via_no_facl", sd_via_no_facl)
        write(file, "sd_via_no_selff", sd_via_no_selff)
        write(file, "sd_via_no_selfc", sd_via_no_selfc)
        write(file, "sd_st_comp", sd_st_comp)
        write(file, "sd_st_facl", sd_st_facl)
        write(file, "sd_st_selfc", sd_st_selfc)
        write(file, "sd_st_selff", sd_st_selff)
    end
    return (nothing)
end

@time intstats()
