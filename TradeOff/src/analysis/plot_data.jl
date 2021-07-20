# A script to read in and analyse model output data
using TradeOff
using JLD
using Plots
import PyPlot

# function to plot the trajectories
function plot_traj()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 2
        error("insufficent inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    rN = 0
    ims = 0
    # Check that all arguments can be converted to integers
    try
        rN = parse(Int64,ARGS[1])
        ims = parse(Int64,ARGS[2])
    catch e
            error("need to provide 2 integers")
    end
    println("Compiled")
    # Extract other simulation parameters from the function
    Np, Rls, Rus, Nt, M = sim_paras()
    # Read in appropriate files
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Species/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    ofile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Species/Run$(rN)Data$(ims)Ims.jld"
    if ~isfile(ofile)
        error("$(ims) immigrations run $(rN) is missing an output file")
    end
    # Read in relevant data
    ps = load(pfile,"ps")
    traj = load(ofile,"traj")
    T = load(ofile,"T")
    micd = load(ofile,"micd")
    its = load(ofile,"its")
    println("Data read in")
    # Find C from a function
    C = merge_data(ps,traj,T,micd,its)
    println("Data merged")
    # Time snapshots that I'm highlighting
    times = [5e6,1.5e7,T[end]]
    # Find total number of strains
    totN = length(micd)
    pyplot(dpi=200)
    # Plot all the populations
    p1 = plot(yaxis=:log10,ylabel="Population (# cells)",ylims=(1e-5,Inf))
    for i = 1:totN
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0)
        plot!(p1,T[inds],C[inds,i],label="")
    end
    vline!(p1,times,color=:red,style=:dash,label="")
    savefig(p1,"Output/pops.png")
    # Plot all the concentrations
    p2 = plot(yaxis=:log10,ylabel="Concentration")#,ylims=(1e-15,Inf))
    for i = 1:ps.M
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,totN+i] .> 0)
        plot!(p2,T[inds],C[inds,totN+i],label="")
    end
    savefig(p2,"Output/concs.png")
    plot(T,C[:,(totN+ps.M+1):(2*totN+ps.M)],label="")
    savefig("Output/as.png")
    plot(T,C[:,(2*totN+ps.M+1):end],label="")
    savefig("Output/fracs.png")
    return(nothing)
end

# Function to plot the generalist vs specialist trade-off with time
function plot_genvsspec()
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
    # Load in hardcoded simulation parameters
    Np, Rls, Rus, Nt, M = sim_paras()
    # Read in parameter file
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Species/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    # Load parameters
    ps = load(pfile,"ps")
    # Define list of times to check
    times = [0.0,5e6,1e7,1.5e7,2e7]
    # Preallocate containers for the reaction distributions
    Rds = zeros(5,length(times)+1)
    # And corresponding biomass distributions
    Rbs = zeros(5,length(times)+1)
    # List of pools already loaded in
    pls = []
    # Array of array to store pools
    pools = Array{Array{Microbe,1},1}(undef,1)
    # Loop over number of repeats
    for i = 1:rps
        # Load in relevant output file
        ofile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Species/Run$(i)Data$(ims)Ims.jld"
        if ~isfile(ofile)
            error("$(ims) immigrations run $(rN) is missing an output file")
        end
        # Load in microbe data, and immigration times
        T = load(ofile,"T")
        traj = load(ofile,"traj")
        micd = load(ofile,"micd")
        its = load(ofile,"its")
        # Find final time
        Tf = T[end]
        # Check if final time is below for any run
        if times[end] > Tf
            println("Run $(i) final time = $(Tf)")
        end
        # Make vector of all times to check
        Ts = cat(times,Tf,dims=1)
        # Now loop over all times to check
        for j = 1:length(Ts)
            # Find initial microbes
            iis = ((micd.↦:ImT) .<= Ts[j]) .& (((micd.↦:ExT) .> Ts[j]) .| isnan.(micd.↦:ExT))
            # Extract subset of the data
            imicd = micd[iis]
            # Find immigration number that we are interested in
            imN = findfirst(x->x>=Ts[j],its)
            # Now extract trajectories we are interested in
            if !isnothing(imN)
                trc = traj[imN]
            else
                # Find current size
                sz = (size(traj[end],2) - ps.M)/3
                sz = convert(Int64,sz)
                # find indices that don't end in zero
                inds = (traj[end])[end,1:sz] .> 0.0
                # Make into the require index form for full traj object
                tinds = collect(1:sz)[inds]
                # Extract relvant trajectory here
                trc = (traj[end])[:,tinds]
            end
            # Then loop over the data
            for k = 1:length(imicd)
                # check for case where pool hasn't already been loaded in
                if imicd[k].PID ∉ pls
                    # Add new pool ID in
                    pls = cat(pls,imicd[k].PID,dims=1)
                    # Find name of pool
                    file = "Pools/ID=$(imicd[k].PID)N=$(Nt)M=$(ps.M)Reacs$(Rls[1])-$(Rus[1]).jld"
                    # Check if this is the first pool
                    if length(pls) == 1
                        # If so save the pool
                        pools[1] = load(file,"mics")
                    else
                        # Otherwise just cat it on existing vector
                        pool = load(file,"mics")
                        pools = cat(pools,pool,dims=1)
                    end
                end
                # Find correct pool to read from
                ind = findfirst(x->x==imicd[k].PID,pls)
                # Use this index to find the correct microbe
                mic = (pools[ind])[imicd[k].MID]
                # Add one to the relevant total
                Rds[mic.R,j] += 1
                # Find biomass of this microbe at end of this time
                if j != 1
                    bm = trc[end,k]
                else
                    # Want starting amounts for inital values
                    bm = trc[1,k]
                end
                # Add biomass to the relevant biomass totals
                Rbs[mic.R,j] += bm
            end
        end
    end
    # Extract maximum reaction number from the pool
    maxR, _ = findmax(pools[1].↦:R)
    # Preallocate output
    pRds = zeros(maxR)
    # Loop of reaction number
    for i = 1:maxR
        # Find all microbes with i reactions
        inds = findall(x->x==i,pools[1].↦:R)
        # Then store the number
        pRds[i] = length(inds)
    end
    pyplot(dpi=200)
    # First plot the distribution in the pool
    plot(ylabel="Number of strains",xlabel="Number of reactions",title="Original pool")
    bar!(pRds,label="")
    savefig("Output/ReacsInitialPool.png")
    for i = 1:length(times)
        plot(ylabel="Number of strains",xlabel="Number of reactions",title="Ecosystem at time $(times[i])s")
        bar!(Rds[:,i],label="")
        savefig("Output/ReacsT=$(i).png")
    end
    # Save final time
    plot(ylabel="Number of strains",xlabel="Number of reactions",title="Final ecosystem")
    bar!(Rds[:,end],label="")
    savefig("Output/ReacsT=$(length(times)+1).png")
    # Same plots for biomass
    for i = 1:length(times)
        plot(ylabel="Strain biomass",xlabel="Number of reactions",title="Ecosystem at time $(times[i])s")
        bar!(Rbs[:,i],label="")
        savefig("Output/BiomassT=$(i).png")
    end
    # Save final time
    plot(ylabel="Strain biomass",xlabel="Number of reactions",title="Final ecosystem")
    bar!(Rbs[:,end],label="")
    savefig("Output/BiomassT=$(length(times)+1).png")
    # Maybe also want average biomass per strain
    for i = 1:length(times)
        plot(ylabel="Average biomass per strain",xlabel="Number of reactions",title="Ecosystem at time $(times[i])s")
        bar!(Rbs[:,i]./Rds[:,i],label="")
        savefig("Output/AvBiomT=$(i).png")
    end
    # Save final time
    plot(ylabel="Average biomass per strain",xlabel="Number of reactions",title="Final ecosystem")
    bar!(Rbs[:,end]./Rds[:,end],label="")
    savefig("Output/AvBiomT=$(length(times)+1).png")
    return(nothing)
end

# Function to load in and plot the averages
function plot_aves()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 1
        error("insufficent inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    ims = 0
    # Check that all arguments can be converted to integers
    try
        ims = parse(Int64,ARGS[1])
    catch e
            error("need to provide an integer")
    end
    println("Compiled")
    # Load in hardcoded simulation parameters
    Np, Rls, Rus, Nt, M = sim_paras()
    # Find file name to load in
    sfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Species/RunStats$(ims)Ims.jld"
    # Check it actually exists
    if ~isfile(sfile)
        error("missing stats file for $(ims) immigrations simulations")
    end
    # Now load out the times, and number of trajectories
    times = load(sfile,"times")
    no_sims = load(sfile,"no_sims")
    # Load in averages
    mn_svt = load(sfile,"mn_svt")
    mn_tsvt = load(sfile,"mn_tsvt")
    mn_sbs = load(sfile,"mn_sbs")
    mn_Rs = load(sfile,"mn_Rs")
    mn_via_R = load(sfile,"mn_via_R")
    mn_ηs_R = load(sfile,"mn_ηs_R")
    mn_ηs = load(sfile,"mn_ηs")
    mn_via_η = load(sfile,"mn_via_η")
    mn_fr_ΔG = load(sfile,"mn_fr_ΔG")
    # Load in standard deviations
    sd_svt = load(sfile,"sd_svt")
    sd_tsvt = load(sfile,"sd_tsvt")
    sd_sbs = load(sfile,"sd_sbs")
    sd_Rs = load(sfile,"sd_Rs")
    sd_via_R = load(sfile,"sd_via_R")
    sd_ηs_R = load(sfile,"sd_ηs_R")
    sd_ηs = load(sfile,"sd_ηs")
    sd_via_η = load(sfile,"sd_via_η")
    sd_fr_ΔG = load(sfile,"sd_fr_ΔG")
    # Preallocate standard errors
    se_Rs = zeros(size(sd_Rs))
    se_via_R = zeros(size(sd_via_R))
    se_ηs_R = zeros(size(sd_ηs_R))
    # Calculate standard errors from this
    se_svt = sd_svt./sqrt.(no_sims)
    se_tsvt = sd_tsvt./sqrt.(no_sims)
    se_sbs = sd_sbs./sqrt.(no_sims)
    se_ηs = sd_ηs./sqrt.(no_sims)
    se_via_η = sd_via_η./sqrt.(no_sims)
    se_fr_ΔG = sd_fr_ΔG./sqrt.(no_sims)
    for i = 1:size(sd_Rs,1)
        se_Rs[i,:] = sd_Rs[i,:]./sqrt.(no_sims)
        se_via_R[i,:] = sd_via_R[i,:]./sqrt.(no_sims)
        se_ηs_R[i,:] = sd_ηs_R[i,:]./sqrt.(no_sims)
    end
    # Setup plotting
    pyplot(dpi=200)
    plot(xlabel="Time (s)",ylabel="Number of surviving strains",xlim=(-Inf,5e7))
    plot!(times,mn_svt,ribbon=se_svt,label="")
    savefig("Output/AvSurvTime.png")
    plot(xlabel="Time (s)",ylabel="Number of viable strains",xlim=(-Inf,5e7))
    plot!(times,mn_tsvt,ribbon=se_tsvt,label="")
    savefig("Output/AvViaTime.png")
    plot(xlabel="Time (s)",ylabel="Number of diversified substrates",xlim=(-Inf,5e7))
    plot!(times,mn_sbs,ribbon=se_sbs,label="")
    savefig("Output/AvSubTime.png")
    plot(xlabel="Time (s)",ylabel="Number of strains",xlim=(-Inf,5e7))
    plot!(times,mn_Rs[1,:],ribbon=se_Rs[1,:],label="R=1")
    plot!(times,mn_Rs[2,:],ribbon=se_Rs[2,:],label="R=2")
    plot!(times,mn_Rs[3,:],ribbon=se_Rs[3,:],label="R=3")
    plot!(times,mn_Rs[4,:],ribbon=se_Rs[4,:],label="R=4")
    plot!(times,mn_Rs[5,:],ribbon=se_Rs[5,:],label="R=5")
    savefig("Output/AvReacsTime.png")
    plot(xlabel="Time (s)",ylabel="Number of strains",xlim=(-Inf,5e7))
    plot!(times,mn_via_R[1,:],ribbon=se_via_R[1,:],label="R=1")
    plot!(times,mn_via_R[2,:],ribbon=se_via_R[2,:],label="R=2")
    plot!(times,mn_via_R[3,:],ribbon=se_via_R[3,:],label="R=3")
    plot!(times,mn_via_R[4,:],ribbon=se_via_R[4,:],label="R=4")
    plot!(times,mn_via_R[5,:],ribbon=se_via_R[5,:],label="R=5")
    savefig("Output/AvViaReacsTime.png")
    plot(xlabel="Time (s)",ylabel="Average eta value",xlim=(-Inf,5e7))
    plot!(times,mn_ηs_R[1,:],ribbon=se_ηs_R[1,:],label="R=1")
    plot!(times,mn_ηs_R[2,:],ribbon=se_ηs_R[2,:],label="R=2")
    plot!(times,mn_ηs_R[3,:],ribbon=se_ηs_R[3,:],label="R=3")
    plot!(times,mn_ηs_R[4,:],ribbon=se_ηs_R[4,:],label="R=4")
    plot!(times,mn_ηs_R[5,:],ribbon=se_ηs_R[5,:],label="R=5")
    savefig("Output/AvEtaperReacTime.png")
    plot(xlabel="Time (s)",ylabel="Average eta value",xlim=(-Inf,5e7))
    plot!(times,mn_ηs,ribbon=se_ηs,label="")
    savefig("Output/AvEtaTime.png")
    plot(xlabel="Time (s)",ylabel="Average eta value",xlim=(-Inf,5e7))
    plot!(times,mn_via_η,ribbon=se_via_η,label="")
    savefig("Output/AvViaEtaTime.png")
    plot(xlabel="Time (s)",ylabel="Fraction transduced",xlim=(-Inf,5e7))
    plot!(times,mn_fr_ΔG,ribbon=se_fr_ΔG,label="")
    savefig("Output/AvFracTransTime.png")
    return(nothing)
end

# function to plot the averages for just one run
function plot_run_averages()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 2
        error("insufficent inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    ims = 0
    rN = 0
    # Check that all arguments can be converted to integers
    try
        ims = parse(Int64,ARGS[1])
        rN = parse(Int64,ARGS[2])
    catch e
            error("need to provide two integers")
    end
    println("Compiled")
    # Load in hardcoded simulation parameters
    Np, Rls, Rus, Nt, M = sim_paras()
    # Read in appropriate files
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Species/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    ofile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Species/AvRun$(rN)Data$(ims)Ims.jld"
    if ~isfile(ofile)
        error("$(ims) immigrations run $(rN) is missing an output file")
    end
    # Read in relevant data
    ps = load(pfile,"ps")
    T = load(ofile,"T")
    svt = load(ofile,"svt")
    ηs = load(ofile,"ηs")
    no_comp = load(ofile,"no_comp")
    no_facl = load(ofile,"no_facl")
    no_self = load(ofile,"no_self",)
    st_comp = load(ofile,"st_comp")
    st_facl = load(ofile,"st_facl")
    st_self = load(ofile,"st_self")
    # Setup plotting
    pyplot(dpi=200)
    # Plot this data
    plot(T,svt,label="",xlabel="Time (s)",ylabel="Number of survivors")
    savefig("Output/SvTime.png")
    plot(T,ηs,label="",xlabel="Time (s)",ylabel="eta")
    savefig("Output/EtaTime.png")
    plot(T,no_comp,label="",xlabel="Time (s)",ylabel="Number of competitive interactions")
    savefig("Output/NoCompTime.png")
    plot(T,no_facl,label="",xlabel="Time (s)",ylabel="Number of facilitation interactions")
    savefig("Output/NoFaclTime.png")
    plot(T,no_self,label="",xlabel="Time (s)",ylabel="Number of self-facilitation interactions")
    savefig("Output/NoSelfTime.png")
    plot(T,log10.(st_comp),label="",xlabel="Time (s)",ylabel="Strength of competitive interactions")
    savefig("Output/StCompTime.png")
    plot(T,log10.(st_facl),label="",xlabel="Time (s)",ylabel="Strength of facilitation interactions")
    savefig("Output/StFaclTime.png")
    plot(T,log10.(st_self),label="",xlabel="Time (s)",ylabel="Strength of self-facilitation interactions")
    savefig("Output/StSelfTime.png")
    return(nothing)
end

# @time plot_run_averages()
@time plot_aves()
# @time plot_traj()
# @time plot_genvsspec()
