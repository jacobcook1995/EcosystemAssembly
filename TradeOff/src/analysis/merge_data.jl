# Script to allow merging of data into the form of full plotable trajectories
export merge_data, sim_paras

# Hard code simulation parameters into this function
function sim_paras()
    # Set the hardcoded variables here
    Np = 1
    Rls = [1]
    Rus = [5]
    Nt = 5000
    M = 25
    return(Np,Rls,Rus,Nt,M)
end

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
