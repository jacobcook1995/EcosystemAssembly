# Script to find how variables change over time, which then saves them
using TradeOff
using JLD

# Function to calculate Shannon diversity from a vector of populations
function shan(pops::Array{Float64, 1})
    # Set initial value
    H = 0.0
    # Make vector of relative abundances
    prbs = pops ./ sum(pops)
    # Loop over relative abundances
    for i in eachindex(prbs)
        # Calculate Shannon entropy for each point
        H -= prbs[i] * log(prbs[i])
    end
    return (H)
end

function v_over_t()
    # Check that sufficient arguments have been provided
    if length(ARGS) < 3
        error("insufficient inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    rps = 0
    ims = 0
    sim_type = 0
    # Check that all arguments can be converted to integers
    try
        rps = parse(Int64, ARGS[1])
        ims = parse(Int64, ARGS[2])
        sim_type = parse(Int64, ARGS[3])
    catch e
        error("need to provide 3 integers")
    end
    println("Compiled")
    flush(stdout)
    # Load in hardcoded simulation parameters
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Token to insert into filenames
    tk = ""
    # Overwritten for no immigration case
    if sim_type == 5
        tk = "NoImm"
    end
    # Read in parameter file
    pfile = "Output/$(tk)$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Paras$(ims)Ims.jld"
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
        ofile = "Output/$(tk)$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Run$(i)Data$(ims)Ims.jld"
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
                file = "Pools/ID=$(micd[j].PID)N=$(Nt)M=$(ps.M)d=$(d)u=$(μrange).jld"
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
        # Preallocate containers to store number of survivors with time
        surviving_species = Array{Int64, 1}(undef, length(T))
        viable_species = Array{Int64, 1}(undef, length(T))
        total_population = zeros(length(T))
        shannon_diversity = zeros(length(T))
        no_substrates = Array{Int64, 1}(undef, length(T))
        species_per_reac_class = Array{Int64, 2}(undef, NoR, length(T))
        viable_species_per_reac_class = Array{Int64, 2}(undef, NoR, length(T))
        average_no_reac_steps = zeros(length(T))
        average_η = zeros(length(T))
        average_ω = zeros(length(T))
        total_biomass_of_viable_species = zeros(length(T))
        average_ΔG = zeros(length(T))
        average_η_per_reac_class = zeros(NoR, length(T))
        average_KS_per_reac_class = zeros(NoR, length(T))
        # Save total number of strains
        total_species = length(micd)
        # Make vector of indices
        ϕ_i = collect((2 * total_species + ps.M + 1):(3 * total_species + ps.M))
        # Loop over all time points
        for j in 1:length(T)
            # Find indices of surviving strains
            inds = findall(x -> x > 1e-5, C[j, 1:total_species])
            # Save number of surviving strains at each time point
            surviving_species[j] = length(inds)
            # Check if at least one species survives
            if surviving_species[j] > 0
                # Save the total population at this point
                total_population[j] = sum(C[j, inds])
                # Find diversity of this point
                shannon_diversity[j] = shan(C[j, inds])
            end
            # Find indices of "viable" species
            vinds = findall(x -> x > 1e5, C[j, 1:total_species])
            # Save number of "viable" species
            viable_species[j] = length(vinds)
            # Then also number of substrates (skipping final waste product)
            no_substrates[j] = count(x -> x > 1e-12,
                                     C[j, (total_species + 1):(total_species + ps.M - 1)])
            # Loop over number of reactions
            for k in 1:NoR
                # Count number of strains with reaction for each case
                species_per_reac_class[k, j] = count(x -> x == k, ms[inds] .↦ :R)
                viable_species_per_reac_class[k, j] = count(x -> x == k,
                                                            ms[vinds] .↦ :R)
            end
            # Set up counters for the number of species with each reaction gap
            c = zeros(M - 1)
            # Find (weighted) total eta value for viable species
            for k in eachindex(vinds)
                average_η[j] += sum(ms[vinds[k]].η .* ms[vinds[k]].ϕP) * C[j, vinds[k]]
                average_ω[j] += ms[vinds[k]].ω * C[j, vinds[k]]
                total_biomass_of_viable_species[j] += C[j, vinds[k]]
                # Loop over all reactions this strain has
                for l in 1:(ms[vinds[k]].R)
                    # Find reaction number
                    Rn = ms[vinds[k]].Reacs[l]
                    # Find relevant reaction
                    r = ps.reacs[Rn]
                    # Then calculate frac transduced
                    average_ΔG[j] += C[j, vinds[k]] * ms[vinds[k]].η[l] .*
                                     ms[vinds[k]].ϕP[l] * ΔGATP / (-r.ΔG0)
                    # Find step size
                    s_size = (r.Prd - r.Rct)
                    # weight this step size to 1 and add to total
                    average_no_reac_steps[j] += C[j, vinds[k]] * (s_size) *
                                                ms[vinds[k]].ϕP[l]
                end
            end
            # Average over biomass of viable species
            if viable_species[j] > 0
                average_η[j] /= total_biomass_of_viable_species[j]
                average_ω[j] /= total_biomass_of_viable_species[j]
                average_no_reac_steps[j] /= total_biomass_of_viable_species[j]
                average_ΔG[j] /= total_biomass_of_viable_species[j]
            end
            # Break down eta and omega value by R
            for k in eachindex(vinds)
                # Find relevant reaction number
                l = ms[vinds[k]].R
                # Add contribution to relevant total
                average_η_per_reac_class[l, j] += sum(ms[vinds[k]].η .* ms[vinds[k]].ϕP)
                average_KS_per_reac_class[l, j] += sum(ms[vinds[k]].KS .* ms[vinds[k]].ϕP)
            end
            # Now weight by number of strains with each type of reaction
            for k in 1:NoR
                if viable_species_per_reac_class[k, j] > 0
                    average_η_per_reac_class[k, j] /= viable_species_per_reac_class[k, j]
                    average_KS_per_reac_class[k, j] /= viable_species_per_reac_class[k, j]
                end
            end
        end
        # Preallocate final ϕR values
        final_ϕR = zeros(surviving_species[end])
        # Find indices of all surviving species
        inds = findall(x -> x > 1e-5, C[end, 1:total_species])
        # Loop over number of survivors
        for j in 1:surviving_species[end]
            # Store corresponding final ϕR value
            final_ϕR[j] = C[end, ϕ_i[inds[j]]]
        end
        # Now just save the relevant data
        jldopen("Output/$(tk)$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/AvRun$(i)Data$(ims)Ims.jld",
                "w") do file
            # Save full time course
            write(file, "T", T)
            # Save reaction data
            write(file, "species_per_reac_class", species_per_reac_class)
            write(file, "viable_species_per_reac_class", viable_species_per_reac_class)
            write(file, "average_η_per_reac_class", average_η_per_reac_class)
            write(file, "average_KS_per_reac_class", average_KS_per_reac_class)
            # Save the other quantities
            write(file, "surviving_species", surviving_species)
            write(file, "viable_species", viable_species)
            write(file, "total_population", total_population)
            write(file, "total_biomass_of_viable_species", total_biomass_of_viable_species)
            write(file, "shannon_diversity", shannon_diversity)
            write(file, "no_substrates", no_substrates)
            write(file, "average_η", average_η)
            write(file, "average_no_reac_steps", average_no_reac_steps)
            write(file, "average_ω", average_ω)
            write(file, "average_ΔG", average_ΔG)
            write(file, "final_ϕR", final_ϕR)
            # Finally save final time to help with benchmarking
            write(file, "final_time_point", T[end])
        end
        println("Run $i analysed")
        flush(stdout)
    end
    return (nothing)
end

@time v_over_t()
