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
    # Find total number of strains
    totN = length(micd)
    pyplot(dpi=200)
    # Plot all the populations
    p1 = plot(yaxis=:log10,ylabel="Population (# cells)")
    for i = 1:totN
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0)
        plot!(p1,T[inds],C[inds,i],label="")
    end
    savefig(p1,"Output/pops.png")
    plot(T,C[:,(totN+1):(totN+ps.M)],label="")
    savefig("Output/concs.png")
    plot(T,C[:,(totN+ps.M+1):(2*totN+ps.M)],label="")
    savefig("Output/as.png")
    plot(T,C[:,(2*totN+ps.M+1):end],label="")
    savefig("Output/fracs.png")
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
