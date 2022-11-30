# Script to find how variables change over time, which then saves them
using TradeOff
using JLD

# function to calculate the strength of facilitation through a particular metabolite
function facl(Mn::Int64,
              ps::TOParameters,
              ms::Array{Microbe, 1},
              pops::Array{Float64, 1},
              concs::Array{Float64, 1},
              ϕRs::Array{Float64, 1},
              fcls::Array{Int64, 2})
    # Find all reactions where metabolite is a product
    Pind = findall(x -> x == Mn, ps.reacs .↦ :Prd)
    # Preallocate counters for self-facilitation and facilitation
    selff = 0
    fcl = 0
    # Now loop over all pairs of strains (strain i facilitates strain j)
    for i in 1:length(ms)
        for j in 1:length(ms)
            # Check that facilitation interaction is present
            if fcls[i, j] > 0
                if i == j
                    # Loop over all of strain i's reaction
                    for k in 1:(ms[i].R)
                        # Only use the ones that count
                        if ms[i].Reacs[k] ∈ Pind
                            # Save reaction
                            r = ps.reacs[ms[i].Reacs[k]]
                            # Find enzyme dedicated to the reaction
                            E = Eα(ϕRs[i], ms[i], k)
                            # Finally find flux for this reaction
                            selff += pops[i] *
                                     qs(concs[r.Rct], concs[r.Prd], E, k, ms[i], ps.T, r) /
                                     NA
                        end
                    end
                else
                    # Loop over all of strain i's reaction
                    for k in 1:(ms[i].R)
                        # Only use the ones that count
                        if ms[i].Reacs[k] ∈ Pind
                            # Save reaction
                            r = ps.reacs[ms[i].Reacs[k]]
                            # Find enzyme dedicated to the reaction
                            E = Eα(ϕRs[i], ms[i], k)
                            # Finally find flux for this reaction
                            fcl += pops[i] *
                                   qs(concs[r.Rct], concs[r.Prd], E, k, ms[i], ps.T, r) / NA
                        end
                    end
                end
            end
        end
    end
    return (selff, fcl)
end

# function to calculate the strength of competition via a particular metabolite
function comp(Mn::Int64,
              ps::TOParameters,
              ms::Array{Microbe, 1},
              pops::Array{Float64, 1},
              concs::Array{Float64, 1},
              ϕRs::Array{Float64, 1},
              cmps::Array{Int64, 2})
    # Find all reactions where metabolite is a substrate
    Sind = findall(x -> x == Mn, ps.reacs .↦ :Rct)
    # Preallocate counters for self-facilitation and facilitation
    selfc = 0
    cmp = 0
    # Now loop over all pairs of strains (strain i facilitates strain j)
    # Strain i's competition with strain j is set by how much strain i consumes
    for i in 1:length(ms)
        for j in 1:length(ms)
            # Find self interactions
            if i == j
                # Loop over all of strain i's reaction
                for k in 1:(ms[i].R)
                    # Save reaction
                    r = ps.reacs[ms[i].Reacs[k]]
                    # Find enzyme dedicated to the reaction
                    E = Eα(ϕRs[i], ms[i], k)
                    # Finally find flux for this reaction
                    selfc += pops[i] *
                             qs(concs[r.Rct], concs[r.Prd], E, k, ms[i], ps.T, r) / NA
                end
                # Look for competition between strains now
            elseif cmps[i, j] > 0
                # Loop over all of strain i's reaction
                for k in 1:(ms[i].R)
                    # Only use the ones that count
                    if ms[i].Reacs[k] ∈ Sind
                        # Save reaction
                        r = ps.reacs[ms[i].Reacs[k]]
                        # Find enzyme dedicated to the reaction
                        E = Eα(ϕRs[i], ms[i], k)
                        # Finally find flux for this reaction
                        cmp += pops[i] *
                               qs(concs[r.Rct], concs[r.Prd], E, k, ms[i], ps.T, r) /
                               NA
                    end
                end
            end
        end
    end
    return (selfc, cmp)
end

# Calculate interaction strengths over time
function ints_over_t()
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
    # Read in parameter file
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations is missing a parameter file")
    end
    # Load parameters
    ps = load(pfile, "ps")
    # List of pools already loaded in
    pls = []
    # Array of array to store pools
    pools = Array{Array{Microbe, 1}, 1}(undef, 1)
    # Counter for number of reactions
    NoR = 0
    # Loop over number of repeats
    for i in 1:rps
        # Load in relevant output file
        ofile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/Run$(i)Data$(ims)Ims.jld"
        if ~isfile(ofile)
            error("$(ims) immigrations run $(rN) is missing an output file")
        end
        # Load in microbe data, and immigration times
        T = load(ofile, "T")
        traj = load(ofile, "traj")
        micd = load(ofile, "micd")
        its = load(ofile, "its")
        # Use to construct full trajectory C
        C = merge_data(ps, traj, T, micd, its)
        # Preallocate vector of microbes
        ms = Array{Microbe, 1}(undef, length(micd))
        # Loop over and find each one
        for j in eachindex(micd)
            # check for case where pool hasn't already been loaded in
            if micd[j].PID ∉ pls
                # Add new pool ID in
                pls = cat(pls, micd[j].PID, dims = 1)
                # Find name of pool
                file = "Pools/ID=$(micd[j].PID)N=$(Nt)M=$(ps.M)d=$(d).jld"
                # Check if this is the first pool
                if length(pls) == 1
                    # If so save the pool
                    pools[1] = load(file, "mics")
                    # Find number of reactions based on this
                    NoR = maximum(pools[1] .↦ :R)
                else
                    # Otherwise just cat it on existing vector
                    pool = load(file, "mics")
                    pools = cat(pools, pool, dims = 1)
                    # Find maximum number of reactions for this pool
                    NoRt = maximum(pools[1] .↦ :R)
                    # Save if higher than old number of reactions
                    NoR = max(NoR, NoRt)
                end
            end
            # Find correct pool to read from
            ind = findfirst(x -> x == micd[j].PID, pls)
            # Use this index to find and save the correct microbe
            ms[j] = (pools[ind])[micd[j].MID]
        end
        # Preallocate interaction matrices
        cmps = zeros(Int64, length(ms), length(ms))
        fcls = zeros(Int64, length(ms), length(ms))
        # Loop over all microbes to make these interaction structure matrices
        for j in eachindex(ms)
            # Loop over microbes to find facilitation terms
            for k in eachindex(ms)
                # Loop over reactions for both strains
                for l in 1:(ms[j].R)
                    for m in 1:(ms[k].R)
                        # Check for facilitation cases
                        if ps.reacs[ms[j].Reacs[l]].Prd == ps.reacs[ms[k].Reacs[m]].Rct
                            fcls[j, k] += 1
                        end
                        # Do the same check for competition cases (avoiding self terms)
                        if j != k &&
                           ps.reacs[ms[j].Reacs[l]].Rct == ps.reacs[ms[k].Reacs[m]].Rct
                            cmps[j, k] += 1
                        end
                    end
                end
            end
        end
        # Preallocate containers to store data of interest with time
        svt = Array{Int64, 1}(undef, length(T))
        tsvt = Array{Int64, 1}(undef, length(T))
        no_comp = Array{Int64, 1}(undef, length(T))
        no_facl = Array{Int64, 1}(undef, length(T))
        no_selff = zeros(Int64, length(T))
        no_selfc = zeros(Int64, length(T))
        via_no_comp = Array{Int64, 1}(undef, length(T))
        via_no_facl = Array{Int64, 1}(undef, length(T))
        via_no_selff = zeros(Int64, length(T))
        via_no_selfc = zeros(Int64, length(T))
        st_comp = zeros(length(T))
        st_selfc = zeros(length(T))
        st_facl = zeros(length(T))
        st_selff = zeros(length(T))
        # Save total number of strains
        numS = length(micd)
        # Loop over all time points
        for j in 1:length(T)
            # Save all strain populations
            pops = C[j, 1:numS]
            # Save concentrations as well
            concs = C[j, (numS + 1):(numS + ps.M)]
            # Also save ribosome fractions
            ϕRs = C[j, (2 * numS + ps.M + 1):(3 * numS + ps.M)]
            # Find indices of surviving strains
            inds = findall(x -> x > 1e-5, pops)
            # Save number of surviving strains at each time point
            svt[j] = length(inds)
            # Find indices of "viable" strains
            vinds = findall(x -> x > 1e5, pops)
            # Save number of "viable" strains
            tsvt[j] = length(vinds)
            # Interactions find via submatrices of precalculated matrices
            no_comp[j] = sum(cmps[inds, inds])
            # Find self interaction terms
            for k in eachindex(inds)
                no_selff[j] += fcls[inds[k], inds[k]]
            end
            # Find all facilitation terms
            no_facl[j] = sum(fcls[inds, inds])
            # Remove self interactions from this total
            no_facl[j] -= no_selff[j]
            # Find number of self-competition reactions (sum of reaction numbers across strains)
            if length(inds) > 0
                no_selfc[j] = sum(ms[inds] .↦ :R)
            end
            # Do the same but only for viable strains
            via_no_comp[j] = sum(cmps[vinds, vinds])
            # Find self interaction terms
            for k in eachindex(vinds)
                via_no_selff[j] += fcls[vinds[k], vinds[k]]
            end
            # Find all facilitation terms
            via_no_facl[j] = sum(fcls[vinds, vinds])
            # Remove self interactions from this total
            via_no_facl[j] -= via_no_selff[j]
            # Find number of self-competition reactions (sum of reaction numbers across strains)
            if length(vinds) > 0
                via_no_selfc[j] = sum(ms[vinds] .↦ :R)
            end
            # Loop over all metabolites
            for k in 1:(ps.M)
                # Find strength of facilitation interactions (of both types) involving this metabolite
                selff, fcl = facl(k, ps, ms[inds], pops[inds], concs, ϕRs[inds],
                                  fcls[inds, inds])
                st_selff[j] += selff
                st_facl[j] += fcl
                # Find strength of competition interactions (of both types) involving this metabolite
                selfc, cmp = comp(k, ps, ms[inds], pops[inds], concs, ϕRs[inds],
                                  cmps[inds, inds])
                st_selfc[j] += selfc
                st_comp[j] += cmp
            end
        end
        # Now just save the relevant data
        jldopen("Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/IntsRun$(i)Data$(ims)Ims.jld",
                "w") do file
            # Save full time course
            write(file, "T", T)
            # Output number of strains and number of viable strains
            write(file, "svt", svt)
            write(file, "tsvt", tsvt)
            # Save number of interactions
            write(file, "no_comp", no_comp)
            write(file, "no_facl", no_facl)
            write(file, "no_selff", no_selff)
            write(file, "no_selfc", no_selfc)
            write(file, "via_no_comp", via_no_comp)
            write(file, "via_no_facl", via_no_facl)
            write(file, "via_no_selff", via_no_selff)
            write(file, "via_no_selfc", via_no_selfc)
            # Save strengths of interactions
            write(file, "st_comp", st_comp)
            write(file, "st_facl", st_facl)
            write(file, "st_selfc", st_selfc)
            write(file, "st_selff", st_selff)
            # Finally save final time to help with benchmarking
            write(file, "Tf", T[end])
        end
        println("Run $i analysed")
        flush(stdout)
    end
    return (nothing)
end

@time ints_over_t()
