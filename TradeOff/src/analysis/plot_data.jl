# A script to read in and analyse model output data
using TradeOff
using JLD
using Plots
import PyPlot

# MAYBE THIS FUNCTION SHOULD BE MADE GLOBAL AT SOME POINT
# function to merge my output data into a plotable form
function merge_data(ps::TOParameters,traj::Array{Array{Float64,2},1},T::Array{Float64,1},
                    micd::Array{MicData,1},its::Array{Float64,1})
    # Find total number of microbes
    totN = length(micd)
    # Find total number of immigration attempts
    ims = length(its)
    # Preallocate array to store all the trajectory data
    C = Array{Float64,2}(undef,length(T),3*totN+ps.M)
    # Previous immigration time is initally zero
    tp = 0.0
    # Index is 1 to begin with
    ind_tp = 1
    # Then loop over every immigration attempt
    for i = 1:ims
        # Extract relevant trajectory
        tt = traj[i]
        # Find new immigration time
        tn = its[i]
        # Find index of this time in vector T
        ind_tn = ind_tp + size(tt,1) - 2
        # Find strains that exist in this window
        inds = ((micd.↦:ImT) .<= tp) .& (((micd.↦:ExT) .>= tn) .| isnan.(micd.↦:ExT))
        # From this calculate number of strains
        Ns = sum(inds)
        # setup counter
        cnt = 1
        # Firstly save the concentrations
        C[ind_tp:ind_tn,(totN+1):(totN+ps.M)] = tt[1:end-1,(Ns+1):(Ns+ps.M)]
        # loop over total number of strains
        for j = 1:totN
            # Find strains that exist within this window
            if inds[j] == true
                C[ind_tp:ind_tn,j] = tt[1:end-1,cnt]
                C[ind_tp:ind_tn,totN+ps.M+j] = tt[1:end-1,Ns+ps.M+cnt]
                C[ind_tp:ind_tn,2*totN+ps.M+j] = tt[1:end-1,2*Ns+ps.M+cnt]
                # Then increment counter
                cnt += 1
            else
                # Strains that aren't present have their variables set as NaN
                C[ind_tp:ind_tn,j] .= NaN
                C[ind_tp:ind_tn,totN+ps.M+j] .= NaN
                C[ind_tp:ind_tn,2*totN+ps.M+j] .= NaN
            end
        end
        # Finally update previous time
        tp = tn
        # Update previous index to be one higher than final index last time
        ind_tp = ind_tn + 1
    end
    # find strains that go extinct at the very end or survive past it
    inds = (((micd.↦:ExT) .== T[end]) .| isnan.(micd.↦:ExT))
    # Save number of survivors
    Ns = sum(inds)
    # Extract final relevant trajectory
    tt = traj[ims+1]
    # setup counter
    cnt = 1
    # Firstly save the concentrations
    C[ind_tp:end,(totN+1):(totN+ps.M)] = tt[1:end,(Ns+1):(Ns+ps.M)]
    # loop over total number of strains
    for j = 1:totN
        # Find strains that exist within this window
        if inds[j] == true
            C[ind_tp:end,j] = tt[1:end,cnt]
            C[ind_tp:end,totN+ps.M+j] = tt[1:end,Ns+ps.M+cnt]
            C[ind_tp:end,2*totN+ps.M+j] = tt[1:end,2*Ns+ps.M+cnt]
            # Then increment counter
            cnt += 1
        else
            # Strains that aren't present have their variables set as NaN
            C[ind_tp:end,j] .= NaN
            C[ind_tp:end,totN+ps.M+j] .= NaN
            C[ind_tp:end,2*totN+ps.M+j] .= NaN
        end
    end
    return(C)
end

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
    plot(T,C[:,(totN+1):(totN+ps.M)],label="")
    savefig("Output/concs.png")
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
                    file = "Pools/ID=$(imicd[j].PID)N=$(Nt)M=$(ps.M)Reacs$(Rls[1])-$(Rus[1]).jld"
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

@time plot_traj()
# @time plot_genvsspec()
