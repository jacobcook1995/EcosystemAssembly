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
    Np, Nt, M, d = sim_paras()
    # Read in appropriate files
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    ofile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/Run$(rN)Data$(ims)Ims.jld"
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
    # Find total number of strains
    totN = length(micd)
    # Find indices of extinct strains
    ext = isnan.(C[end,1:totN])
    # Invert to find survivors
    svr = .!ext
    pyplot(dpi=200)
    # Plot all the populations
    p1 = plot(yaxis=:log10,ylabel="Population (# cells)",ylims=(1e-5,Inf))
    for i = 1:totN
        if svr[i] == true
            # Find and eliminate zeros so that they can be plotted on a log plot
            inds = (C[:,i] .> 0)
            plot!(p1,T[inds],C[inds,i],label="")
        end
    end
    savefig(p1,"Output/pops.png")
    # Plot all the concentrations
    p2 = plot(yaxis=:log10,ylabel="Concentration")#,ylims=(1e-15,Inf))
    for i = 1:ps.M
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,totN+i] .> 0)
        plot!(p2,T[inds],C[inds,totN+i],label="")
    end
    savefig(p2,"Output/concs.png")
    # Plot energy concentrat
    p3 = plot(ylabel="Energy Concentration")
    for i = 1:totN
        if svr[i] == true
            plot!(p3,T,C[:,totN+ps.M+i],label="")
        end
    end
    savefig(p3,"Output/as.png")
    # Plot all the ribosome fractions
    p4 = plot(ylabel="Ribosome fraction")
    for i = 1:totN
        if svr[i] == true
            plot!(p4,T,C[:,2*totN+ps.M+i],label="")
        end
    end
    savefig(p4,"Output/fracs.png")
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
    Np, Nt, M, d = sim_paras()
    # Find file name to load in
    sfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/RunStats$(ims)Ims.jld"
    # Check it actually exists
    if ~isfile(sfile)
        error("missing stats file for $(ims) immigrations simulations")
    end
    # Now load out the times, and number of trajectories
    times = load(sfile,"times")
    no_sims = load(sfile,"no_sims")
    no_via = load(sfile,"no_via")
    no_rs = load(sfile,"no_rs")
    # Load in averages
    mn_svt = load(sfile,"mn_svt")
    mn_tsvt = load(sfile,"mn_tsvt")
    mn_sbs = load(sfile,"mn_sbs")
    mn_Rs = load(sfile,"mn_Rs")
    mn_via_R = load(sfile,"mn_via_R")
    mn_ηs_R = load(sfile,"mn_ηs_R")
    mn_ωs_R = load(sfile,"mn_ωs_R")
    mn_kc_R = load(sfile,"mn_kc_R")
    mn_KS_R = load(sfile,"mn_KS_R")
    mn_kr_R = load(sfile,"mn_kr_R")
    mn_ηs = load(sfile,"mn_ηs")
    mn_via_η = load(sfile,"mn_via_η")
    mn_ωs = load(sfile,"mn_ωs")
    mn_via_ω = load(sfile,"mn_via_ω")
    mn_fr_ΔG = load(sfile,"mn_fr_ΔG")
    mn_kcs = load(sfile,"mn_kcs")
    mn_KSs = load(sfile,"mn_KSs")
    mn_krs = load(sfile,"mn_krs")
    mn_av_steps = load(sfile,"mn_av_steps")
    # Load in standard deviations
    sd_svt = load(sfile,"sd_svt")
    sd_tsvt = load(sfile,"sd_tsvt")
    sd_sbs = load(sfile,"sd_sbs")
    sd_Rs = load(sfile,"sd_Rs")
    sd_via_R = load(sfile,"sd_via_R")
    sd_ηs_R = load(sfile,"sd_ηs_R")
    sd_ωs_R = load(sfile,"sd_ωs_R")
    sd_kc_R = load(sfile,"sd_kc_R")
    sd_KS_R = load(sfile,"sd_KS_R")
    sd_kr_R = load(sfile,"sd_kr_R")
    sd_ηs = load(sfile,"sd_ηs")
    sd_via_η = load(sfile,"sd_via_η")
    sd_ωs = load(sfile,"sd_ωs")
    sd_via_ω = load(sfile,"sd_via_ω")
    sd_fr_ΔG = load(sfile,"sd_fr_ΔG")
    sd_kcs = load(sfile,"sd_kcs")
    sd_KSs = load(sfile,"sd_KSs")
    sd_krs = load(sfile,"sd_krs")
    sd_av_steps = load(sfile,"sd_av_steps")
    # Preallocate standard errors
    se_Rs = zeros(size(sd_Rs))
    se_via_R = zeros(size(sd_via_R))
    se_ηs_R = zeros(size(sd_ηs_R))
    se_ωs_R = zeros(size(sd_ωs_R))
    se_kc_R = zeros(size(sd_kc_R))
    se_KS_R = zeros(size(sd_KS_R))
    se_kr_R = zeros(size(sd_kr_R))
    # Calculate standard errors from this
    se_svt = sd_svt./sqrt.(no_sims)
    se_tsvt = sd_tsvt./sqrt.(no_sims)
    se_sbs = sd_sbs./sqrt.(no_sims)
    se_ηs = sd_ηs./sqrt.(no_sims)
    se_ωs = sd_ωs./sqrt.(no_sims)
    # Calculation (slightly) different in the viable case
    se_via_η = sd_via_η./sqrt.(no_via)
    se_via_ω = sd_via_ω./sqrt.(no_via)
    se_fr_ΔG = sd_fr_ΔG./sqrt.(no_via)
    se_kcs = sd_kcs./sqrt.(no_via)
    se_KSs = sd_KSs./sqrt.(no_via)
    se_krs = sd_krs./sqrt.(no_via)
    se_av_steps = sd_av_steps./sqrt.(no_via)
    for i = 1:size(sd_Rs,1)
        se_Rs[i,:] = sd_Rs[i,:]./sqrt.(no_sims)
        se_via_R[i,:] = sd_via_R[i,:]./sqrt.(no_via)
        se_ηs_R[i,:] = sd_ηs_R[i,:]./sqrt.(no_rs[i,:])
        se_ωs_R[i,:] = sd_ωs_R[i,:]./sqrt.(no_rs[i,:])
        se_kc_R[i,:] = sd_kc_R[i,:]./sqrt.(no_rs[i,:])
        se_KS_R[i,:] = sd_KS_R[i,:]./sqrt.(no_rs[i,:])
        se_kr_R[i,:] = sd_kr_R[i,:]./sqrt.(no_rs[i,:])
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
    plot!(times,mn_Rs[3,:],ribbon=se_Rs[3,:],label="R=3")
    plot!(times,mn_Rs[5,:],ribbon=se_Rs[5,:],label="R=5")
    plot!(times,mn_Rs[7,:],ribbon=se_Rs[7,:],label="R=7")
    savefig("Output/AvReacsTime.png")
    plot(xlabel="Time (s)",ylabel="Number of strains",xlim=(-Inf,5e7))
    plot!(times,mn_via_R[1,:],ribbon=se_via_R[1,:],label="R=1")
    plot!(times,mn_via_R[3,:],ribbon=se_via_R[3,:],label="R=3")
    plot!(times,mn_via_R[5,:],ribbon=se_via_R[5,:],label="R=5")
    plot!(times,mn_via_R[7,:],ribbon=se_via_R[7,:],label="R=7")
    savefig("Output/AvViaReacsTime.png")
    plot(xlabel="Time (s)",ylabel="Average eta value",xlim=(-Inf,5e7))
    plot!(times,mn_ηs_R[1,:],ribbon=se_ηs_R[1,:],label="R=1")
    plot!(times,mn_ηs_R[3,:],ribbon=se_ηs_R[3,:],label="R=3")
    plot!(times,mn_ηs_R[5,:],ribbon=se_ηs_R[5,:],label="R=5")
    plot!(times,mn_ηs_R[7,:],ribbon=se_ηs_R[7,:],label="R=7")
    savefig("Output/AvEtaperReacTime.png")
    plot(xlabel="Time (s)",ylabel="Average omega value",xlim=(-Inf,5e7))
    plot!(times,mn_ωs_R[1,:],ribbon=se_ωs_R[1,:],label="R=1")
    plot!(times,mn_ωs_R[3,:],ribbon=se_ωs_R[3,:],label="R=3")
    plot!(times,mn_ωs_R[5,:],ribbon=se_ωs_R[5,:],label="R=5")
    plot!(times,mn_ωs_R[7,:],ribbon=se_ωs_R[7,:],label="R=7")
    savefig("Output/AvOmegaperReacTime.png")
    plot(xlabel="Time (s)",ylabel="Average forward rate constant",xlim=(-Inf,5e7))
    plot!(times,mn_kc_R[1,:],ribbon=se_kc_R[1,:],label="R=1")
    plot!(times,mn_kc_R[3,:],ribbon=se_kc_R[3,:],label="R=3")
    plot!(times,mn_kc_R[5,:],ribbon=se_kc_R[5,:],label="R=5")
    plot!(times,mn_kc_R[7,:],ribbon=se_kc_R[7,:],label="R=7")
    savefig("Output/AvkcperReacTime.png")
    plot(xlabel="Time (s)",ylabel="Average half saturation constant",xlim=(-Inf,5e7))
    plot!(times,mn_KS_R[1,:],ribbon=se_KS_R[1,:],label="R=1")
    plot!(times,mn_KS_R[3,:],ribbon=se_KS_R[3,:],label="R=3")
    plot!(times,mn_KS_R[5,:],ribbon=se_KS_R[5,:],label="R=5")
    plot!(times,mn_KS_R[7,:],ribbon=se_KS_R[7,:],label="R=7")
    savefig("Output/AvKSperReacTime.png")
    plot(xlabel="Time (s)",ylabel="Average reversibility factor",xlim=(-Inf,5e7))
    plot!(times,mn_kr_R[1,:],ribbon=se_kr_R[1,:],label="R=1")
    plot!(times,mn_kr_R[3,:],ribbon=se_kr_R[3,:],label="R=3")
    plot!(times,mn_kr_R[5,:],ribbon=se_kr_R[5,:],label="R=5")
    plot!(times,mn_kr_R[7,:],ribbon=se_kr_R[7,:],label="R=7")
    savefig("Output/AvkrperReacTime.png")
    plot(xlabel="Time (s)",ylabel="Average eta value",xlim=(-Inf,5e7))
    plot!(times,mn_ηs,ribbon=se_ηs,label="")
    savefig("Output/AvEtaTime.png")
    plot(xlabel="Time (s)",ylabel="Average eta value",xlim=(-Inf,5e7))
    plot!(times,mn_via_η,ribbon=se_via_η,label="")
    savefig("Output/AvViaEtaTime.png")
    plot(xlabel="Time (s)",ylabel="Average omega value",xlim=(-Inf,5e7))
    plot!(times,mn_ωs,ribbon=se_ωs,label="")
    savefig("Output/AvOmegaTime.png")
    plot(xlabel="Time (s)",ylabel="Average omega value",xlim=(-Inf,5e7))
    plot!(times,mn_via_ω,ribbon=se_via_ω,label="")
    savefig("Output/AvViaOmegaTime.png")
    plot(xlabel="Time (s)",ylabel="Fraction transduced",xlim=(-Inf,5e7))
    plot!(times,mn_fr_ΔG,ribbon=se_fr_ΔG,label="")
    savefig("Output/AvFracTransTime.png")
    plot(xlabel="Time (s)",ylabel="Average forward rate",xlim=(-Inf,5e7))
    plot!(times,mn_kcs,ribbon=se_kcs,label="")
    savefig("Output/AvkcTime.png")
    plot(xlabel="Time (s)",ylabel="Average saturation constant",xlim=(-Inf,5e7))
    plot!(times,mn_KSs,ribbon=se_KSs,label="")
    savefig("Output/AvKSTime.png")
    plot(xlabel="Time (s)",ylabel="Average reversibility factor",xlim=(-Inf,5e7))
    plot!(times,mn_krs,ribbon=se_krs,label="")
    savefig("Output/AvkrTime.png")
    plot(xlabel="Time (s)",ylabel="Average steps in metabolite hierachy",xlim=(-Inf,5e7))
    plot!(times,mn_av_steps,ribbon=se_av_steps,label="")
    savefig("Output/AvStepsTime.png")
    return(nothing)
end

# Function to load in and plot the average interaction strengths
function plot_av_ints()
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
    Np, Nt, M, d = sim_paras()
    # Find file name to load in
    sfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/IntStats$(ims)Ims.jld"
    # Check it actually exists
    if ~isfile(sfile)
        error("missing stats file for $(ims) immigrations simulations")
    end
    # Now load out the times, and number of trajectories
    times = load(sfile,"times")
    no_sims = load(sfile,"no_sims")
    no_via = load(sfile,"no_via")
    # Load in averages
    mn_svt = load(sfile,"mn_svt")
    mn_tsvt = load(sfile,"mn_tsvt")
    mn_no_comp = load(sfile,"mn_no_comp")
    mn_no_facl = load(sfile,"mn_no_facl")
    mn_no_selff = load(sfile,"mn_no_selff")
    mn_no_selfc = load(sfile,"mn_no_selfc")
    mn_via_no_comp = load(sfile,"mn_via_no_comp")
    mn_via_no_facl = load(sfile,"mn_via_no_facl")
    mn_via_no_selff = load(sfile,"mn_via_no_selff")
    mn_via_no_selfc = load(sfile,"mn_via_no_selfc")
    mn_st_comp = load(sfile,"mn_st_comp")
    mn_st_facl = load(sfile,"mn_st_facl")
    mn_st_selfc = load(sfile,"mn_st_selfc")
    mn_st_selff = load(sfile,"mn_st_selff")
    # Load in standard deviations
    sd_svt = load(sfile,"sd_svt")
    sd_tsvt = load(sfile,"sd_tsvt")
    sd_no_comp = load(sfile,"sd_no_comp")
    sd_no_facl = load(sfile,"sd_no_facl")
    sd_no_selff = load(sfile,"sd_no_selff")
    sd_no_selfc = load(sfile,"sd_no_selfc")
    sd_via_no_comp = load(sfile,"sd_via_no_comp")
    sd_via_no_facl = load(sfile,"sd_via_no_facl")
    sd_via_no_selff = load(sfile,"sd_via_no_selff")
    sd_via_no_selfc = load(sfile,"sd_via_no_selfc")
    sd_st_comp = load(sfile,"sd_st_comp")
    sd_st_facl = load(sfile,"sd_st_facl")
    sd_st_selfc = load(sfile,"sd_st_selfc")
    sd_st_selff = load(sfile,"sd_st_selff")
    # Calculate standard errors from this
    se_svt = sd_svt./sqrt.(no_sims)
    se_tsvt = sd_tsvt./sqrt.(no_sims)
    se_no_comp = sd_no_comp./sqrt.(no_sims)
    se_no_facl = sd_no_facl./sqrt.(no_sims)
    se_no_selff = sd_no_selff./sqrt.(no_sims)
    se_no_selfc = sd_no_selfc./sqrt.(no_sims)
    se_st_comp = sd_st_comp./sqrt.(no_sims)
    se_st_facl = sd_st_facl./sqrt.(no_sims)
    se_st_selff = sd_st_selff./sqrt.(no_sims)
    se_st_selfc = sd_st_selfc./sqrt.(no_sims)
    # Calculation (slightly) different in the viable case
    se_via_no_comp = sd_via_no_comp./sqrt.(no_via)
    se_via_no_facl = sd_via_no_facl./sqrt.(no_via)
    se_via_no_selff = sd_via_no_selff./sqrt.(no_via)
    se_via_no_selfc = sd_via_no_selfc./sqrt.(no_via)
    # Setup plotting
    # ALL QUITE MUCKY HERE AT THE MOMENT, CLEAN THESE UP BEFORE I PRESENT THEM
    pyplot(dpi=200)
    plot(xlabel="Time (s)",ylabel="Number of competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_no_comp./(mn_no_selfc.^2),label="")
    savefig("Output/AvCompTime.png")
    plot(xlabel="Time (s)",ylabel="Number of faciliation interactions",xlim=(-Inf,5e7))
    plot!(times,mn_no_facl./(mn_no_selfc.^2),label="")
    savefig("Output/AvFaclTime.png")
    plot(xlabel="Time (s)",ylabel="Number of (self) competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_no_selfc./mn_no_selfc,label="")
    savefig("Output/AvSelfCompTime.png")
    plot(xlabel="Time (s)",ylabel="Number of (self) faciliation interactions",xlim=(-Inf,5e7))
    plot!(times,mn_no_selff./mn_no_selfc,label="")
    savefig("Output/AvSelfFaclTime.png")
    plot(xlabel="Time (s)",ylabel="Number of viable competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_via_no_comp./(mn_via_no_selfc.^2),label="")
    savefig("Output/AvViaCompTime.png")
    plot(xlabel="Time (s)",ylabel="Number of viable faciliation interactions",xlim=(-Inf,5e7))
    plot!(times,mn_via_no_facl./(mn_via_no_selfc.^2),label="")
    savefig("Output/AvViaFaclTime.png")
    plot(xlabel="Time (s)",ylabel="Number of viable (self) competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_via_no_selfc./mn_via_no_selfc,label="")
    savefig("Output/AvViaSelfCompTime.png")
    plot(xlabel="Time (s)",ylabel="Number of viable (self) competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_via_no_selff./mn_via_no_selfc,label="")
    savefig("Output/AvViaSelfFaclTime.png")
    plot(xlabel="Time (s)",ylabel="Strength of competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_st_comp,ribbon=se_st_comp,label="")
    savefig("Output/AvStCompTime.png")
    plot(xlabel="Time (s)",ylabel="Strength of faciliation interactions",xlim=(-Inf,5e7))
    plot!(times,mn_st_facl,ribbon=se_st_facl,label="")
    savefig("Output/AvStFaclTime.png")
    plot(xlabel="Time (s)",ylabel="Strength of (self) competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_st_selfc,ribbon=se_st_selfc,label="")
    savefig("Output/AvStSelfCompTime.png")
    plot(xlabel="Time (s)",ylabel="Strength of (self) competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_st_selff,ribbon=se_st_selff,label="")
    savefig("Output/AvStSelfFaclTime.png")
    plot(xlabel="Time (s)",ylabel="Mean strength of competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_st_comp./mn_no_comp,label="")
    savefig("Output/MnStCompTime.png")
    plot(xlabel="Time (s)",ylabel="Mean strength of faciliation interactions",xlim=(-Inf,5e7))
    plot!(times,mn_st_facl./mn_no_facl,label="")
    savefig("Output/MnStFaclTime.png")
    plot(xlabel="Time (s)",ylabel="Mean strength of (self) competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_st_selfc./mn_no_selfc,label="")
    savefig("Output/MnStSelfCompTime.png")
    plot(xlabel="Time (s)",ylabel="Mean strength of (self) competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_st_selff./mn_no_selff,label="")
    savefig("Output/MnStSelfFaclTime.png")
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
    Np, Nt, M, d = sim_paras()
    # Read in appropriate files
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    ofile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/AvRun$(rN)Data$(ims)Ims.jld"
    if ~isfile(ofile)
        error("$(ims) immigrations run $(rN) is missing an output file")
    end
    # Read in relevant data
    ps = load(pfile,"ps")
    T = load(ofile,"T")
    svt = load(ofile,"svt")
    tsvt = load(ofile,"tsvt")
    ηs = load(ofile,"ηs")
    # Setup plotting
    pyplot(dpi=200)
    # Plot this data
    plot(T,svt,label="",xlabel="Time (s)",ylabel="Number of survivors")
    savefig("Output/SvTime.png")
    plot(T,ηs,label="",xlabel="Time (s)",ylabel="eta")
    savefig("Output/EtaTime.png")
    plot(T,tsvt,label="",xlabel="Time (s)",ylabel="Number of viable strains")
    savefig("Output/ViaTime.png")
    return(nothing)
end

# @time plot_run_averages()
# @time plot_aves()
@time plot_traj()
# @time plot_genvsspec()
# @time plot_av_ints()
