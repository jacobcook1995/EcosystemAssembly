# Script to find how variables change over time, which then saves them
using TradeOff
using JLD

# function to calculate the "interaction strength" between two facilitating strains
function fcl_flx(mic1::Microbe,mic2::Microbe,nI::Int64,concs::Array{Float64,1},md1::Array{Float64,1},
                    md2::Array{Float64,1},ps::TOParameters)
    # Preallocate vector of indices
    inds = zeros(Int64,nI,2)
    # Set up counter
    c = 0
    # Loop over all reactions
    for i = 1:mic1.R
        for j = 1:mic2.R
            if ps.reacs[mic1.Reacs[i]].Prd == ps.reacs[mic2.Reacs[j]].Rct
                # Increment counter
                c += 1
                # Store indicies
                inds[c,1] = i
                inds[c,2] = j
            end
        end
    end
    # Set initial flux to zero
    av_fl = 0.0
    # Loop over number of interactions
    for i = 1:nI
        # Extract reactions
        r1 = ps.reacs[mic1.Reacs[inds[i,1]]]
        r2 = ps.reacs[mic2.Reacs[inds[i,2]]]
        # Extract initial substrate, intermediate product, and final product
        S = concs[r1.Rct]
        P = concs[r1.Prd]
        W = concs[r2.Prd]
        # Calculate enzyme fractions
        E1 = Eα(md1[2],mic1,inds[i,1])
        E2 = Eα(md2[2],mic2,inds[i,2])
        # Calculate net fluxes
        nf1 = md1[1]*qs(S,P,E1,inds[i,1],mic1,ps.T,r1)
        nf2 = md2[1]*qs(P,W,E2,inds[i,2],mic2,ps.T,r2)
        # Sqrt the product and add to the sum
        av_fl += sqrt(nf1*nf2)
    end
    return(av_fl)
end

# function to calculate the "interaction strength" between two competing strains
function fcl_cmp(mic1::Microbe,mic2::Microbe,nI::Int64,concs::Array{Float64,1},md1::Array{Float64,1},
                    md2::Array{Float64,1},ps::TOParameters)
    # Preallocate vector of indices
    inds = zeros(Int64,nI,2)
    # Set up counter
    c = 0
    # Loop over all reactions
    for i = 1:mic1.R
        for j = 1:mic2.R
            if ps.reacs[mic1.Reacs[i]].Rct == ps.reacs[mic2.Reacs[j]].Rct
                # Increment counter
                c += 1
                # Store indicies
                inds[c,1] = i
                inds[c,2] = j
            end
        end
    end
    # Set initial flux to zero
    av_fl = 0.0
    # Loop over number of interactions
    for i = 1:nI
        # Extract reactions
        r1 = ps.reacs[mic1.Reacs[inds[i,1]]]
        r2 = ps.reacs[mic2.Reacs[inds[i,2]]]
        # Extract initial substrate, and the two products (which will often be the same)
        S = concs[r1.Rct]
        P1 = concs[r1.Prd]
        P2 = concs[r2.Prd]
        # Calculate enzyme fractions
        E1 = Eα(md1[2],mic1,inds[i,1])
        E2 = Eα(md2[2],mic2,inds[i,2])
        # Calculate net fluxes
        nf1 = md1[1]*qs(S,P1,E1,inds[i,1],mic1,ps.T,r1)
        nf2 = md2[1]*qs(S,P2,E2,inds[i,2],mic2,ps.T,r2)
        # Sqrt the product and add to the sum
        av_fl += sqrt(nf1*nf2)
    end
    return(av_fl)
end

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
    Np, Rls, Rus, Nt, M = sim_paras()
    # Save number of reactions
    NoR = Rus[1] - Rls[1] + 1
    # Read in parameter file
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Species/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    # Load parameters
    ps = load(pfile,"ps")
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
                file = "Pools/ID=$(micd[j].PID)N=$(Nt)M=$(ps.M)Reacs$(Rls[1])-$(Rus[1]).jld"
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
            ind = findfirst(x->x==micd[j].PID,pls)
            # Use this index to find and save the correct microbe
            ms[j] = (pools[ind])[micd[j].MID]
        end
        # Preallocate interaction matrices
        cmps = zeros(Int64,length(ms),length(ms))
        fcls = zeros(Int64,length(ms),length(ms))
        # Loop over all microbes to make these interaction structure matrices
        for j = 1:length(ms)
            # Loop over microbes to find facilitation terms
            for k = 1:length(ms)
                # Loop over reactions for both strains
                for l = 1:ms[j].R
                    for m = 1:ms[k].R
                        # Check for facilitation cases
                        if ps.reacs[ms[j].Reacs[l]].Prd == ps.reacs[ms[k].Reacs[m]].Rct
                            fcls[j,k] += 1
                        end
                        # Do the same check for competition cases (avoiding double counting)
                        if j < k && ps.reacs[ms[j].Reacs[l]].Rct == ps.reacs[ms[k].Reacs[m]].Rct
                            cmps[j,k] += 1
                        end
                    end
                end
            end
        end
        # Preallocate containers to store number of survivors with time
        svt = Array{Int64,1}(undef,length(T))
        tsvt = Array{Int64,1}(undef,length(T))
        sbs = Array{Int64,1}(undef,length(T))
        Rs = Array{Int64,2}(undef,NoR,length(T))
        ηs = zeros(length(T))
        no_comp = Array{Int64,1}(undef,length(T))
        no_facl = Array{Int64,1}(undef,length(T))
        no_self = zeros(Int64,length(T))
        st_comp = zeros(length(T))
        st_facl = zeros(length(T))
        st_self = zeros(length(T))
        # Save total number of strains
        numS = length(micd)
        # Loop over all time points
        for j = 1:length(T)
            # Find indices of surviving strains
            inds = findall(x->x>1e-5,C[j,1:numS])
            # Save number of surviving strains at each time point
            svt[j] = length(inds)
            # Also top survivors
            tsvt[j] = count(x->x>1e5,C[j,1:numS])
            # Then also number of substrates
            sbs[j] = count(x->x>1e-12,C[j,(numS+1):(numS+ps.M)])
            # Loop over number of reactions
            for k = 1:NoR
                # Count number of strains with reaction for each case
                Rs[k,j] = count(x->x==k,ms[inds].↦:R)
            end
            # Find (weighted) total eta value
            for k = 1:length(inds)
                ηs[j] += sum(ms[inds[k]].η.*ms[inds[k]].ϕP)
            end
            # Average over number of strains
            if svt[j] > 0
                ηs[j] /= svt[j]
            end
            # Interactions find via submatrices of precalulated matrices
            no_comp[j] = sum(cmps[inds,inds])
            # Find self interaction terms
            for k = 1:length(inds)
                no_self[j] += fcls[inds[k],inds[k]]
            end
            # Find all interaction terms
            no_facl[j] = sum(fcls[inds,inds])
            # Remove self interactions from this total
            no_facl[j] -= no_self[j]
            # Now consider strength of self-interactions
            for k = 1:length(inds)
                if fcls[inds[k],inds[k]] != 0
                    # Find relevant microbe data
                    md = [C[j,inds[k]],C[j,ps.M+2*numS+inds[k]]]
                    # Use to calculate contribution to self-facilitation strength
                    st_self[j] += fcl_flx(ms[inds[k]],ms[inds[k]],fcls[inds[k],inds[k]],C[j,(numS+1):(numS+ps.M)],md,md,ps)
                end
            end
            # Same process facilitation interactions in general
            for k = 1:length(inds)
                for l = 1:length(inds)
                    if k != l && fcls[inds[k],inds[l]] != 0
                        # Find relevant microbe data
                        md1 = [C[j,inds[k]],C[j,ps.M+2*numS+inds[k]]]
                        md2 = [C[j,inds[l]],C[j,ps.M+2*numS+inds[l]]]
                        # Use to calculate contribution to facilitation strength
                        st_facl[j] += fcl_flx(ms[inds[k]],ms[inds[l]],fcls[inds[k],inds[l]],C[j,(numS+1):(numS+ps.M)],md1,md2,ps)
                    end
                end
            end
            # Finally the same process for competition
            for k = 1:length(inds)
                for l = (k+1):length(inds)
                    if cmps[inds[k],inds[l]] != 0
                        # Find relevant microbe data
                        md1 = [C[j,inds[k]],C[j,ps.M+2*numS+inds[k]]]
                        md2 = [C[j,inds[l]],C[j,ps.M+2*numS+inds[l]]]
                        # Use to calculate contribution to competition strength
                        st_comp[j] += fcl_cmp(ms[inds[k]],ms[inds[l]],cmps[inds[k],inds[l]],C[j,(numS+1):(numS+ps.M)],md1,md2,ps)
                    end
                end
            end
        end
        # Now just save the relevant data
        jldopen("Output/$(Np)Pools$(M)Metabolites$(Nt)Species/AvRun$(i)Data$(ims)Ims.jld","w") do file
            # Save full timecourse
            write(file,"T",T)
            # Save reaction data
            write(file,"Rs",Rs)
            # Save the other quantities
            write(file,"svt",svt)
            write(file,"tsvt",tsvt)
            write(file,"sbs",sbs)
            write(file,"ηs",ηs)
            write(file,"no_comp",no_comp)
            write(file,"no_facl",no_facl)
            write(file,"no_self",no_self)
            write(file,"st_comp",st_comp)
            write(file,"st_facl",st_facl)
            write(file,"st_self",st_self)
            # Finally save final time to help with benchmarking
            write(file,"Tf",T[end])
        end
        println("Run $i analysed")
        flush(stdout)
    end
    return(nothing)
end


@time v_over_t()
