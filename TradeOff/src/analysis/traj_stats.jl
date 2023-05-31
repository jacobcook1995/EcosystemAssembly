# File to calculate stats for the variables across trajectories
using TradeOff
using JLD

# Function to make a dictionary to store all desired data based on list of names and
# dimensions. This also involves preallocating data where relevant
function make_data_dictonary_from_list(variables_of_interest::Array{
                                                                    Tuple{String, Int64,
                                                                          Bool, Bool}, 1
                                                                    },
                                       repeats::Int64,
                                       no_reactions::Int64,
                                       no_time_points::Int64)
    # Make initially empty dictionary
    data_dict = Dict()
    # Then use a for loop to populate it with names + preallocated memory
    for variable in variables_of_interest
        if variable[2] == 2
            data_dict[variable[1]] = Dict("combined_data" => zeros(repeats, no_time_points),
                                          "run_data" => Float64[], "dims" => 2,
                                          "viable" => variable[3],
                                          "divide_by_R" => variable[4])
        elseif variable[2] == 3
            data_dict[variable[1]] = Dict("combined_data" => zeros(repeats, no_reactions,
                                                                   no_time_points),
                                          "run_data" => Float64[], "dims" => 3,
                                          "viable" => variable[3],
                                          "divide_by_R" => variable[4])
        end
    end

    return data_dict
end

# This function loads in the relevant data and adds it into the data dictionary
function load_trajectory_vars_to_dict!(vfile::String, data_dict::Dict)
    # Loop over every name in the data dictionary
    for variable in keys(data_dict)
        data_dict[variable]["run_data"] = load(vfile, variable)
    end
    return (data_dict)
end

# This functions uses interpolation to take time slices for each variable. This is done to
# ensure that runs can be sensibly compared
function add_to_combined_data!(data_dict::Dict, T::Array{Float64, 1},
                               times::Array{Float64, 1}, Tind::Int64, cnt::Int64, i::Int64)
    # Skip averaging if previous point is missing
    if Tind > 1 && data_dict["viable_species"]["run_data"][Tind - 1] != 0
        # Calculate relevant time gaps
        Tg = (T[Tind] - T[Tind - 1])
        T1x = times[cnt] - T[Tind - 1]
        T2x = T[Tind] - times[cnt]
        # Loop over every name in the data dictionary
        for variable in keys(data_dict)
            # Check variable dimensionality
            if data_dict[variable]["dims"] == 2
                # Interpolate and add into combined data
                point = interpolate_time(data_dict[variable]["run_data"], Tg, T1x, T2x,
                                         Tind)
                data_dict[variable]["combined_data"][i, cnt] = point
            elseif data_dict[variable]["dims"] == 3
                point = interpolate_time(data_dict[variable]["run_data"], Tg, T1x, T2x,
                                         Tind)
                data_dict[variable]["combined_data"][i, :, cnt] = point
            end
        end
    else
        # Loop over every name in the data dictionary
        for variable in keys(data_dict)
            # Check variable dimensionality
            if data_dict[variable]["dims"] == 2
                # In the one case just add the value at time = 0
                point = data_dict[variable]["run_data"][Tind]
                data_dict[variable]["combined_data"][i, cnt] = point
            elseif data_dict[variable]["dims"] == 3
                point = data_dict[variable]["run_data"][:, Tind]
                data_dict[variable]["combined_data"][i, :, cnt] = point
            end
        end
    end
    return (data_dict)
end

# This function calculates the means for all variables along the trajectory
function calculate_trajectory_means!(data_dict::Dict, no_reactions::Int64,
                                     no_simulations::Vector{Float64},
                                     no_viable_simulations::Vector{Float64},
                                     no_simulations_with_R::Array{Float64, 2})
    # Loop over every name in the data dictionary
    for variable in keys(data_dict)
        # Sum to find totals
        total = dropdims(sum(data_dict[variable]["combined_data"], dims = 1), dims = 1)
        # Check variable dimensionality (3D case requires preallocation)
        if data_dict[variable]["dims"] == 2
            # Now calculate means
            if data_dict[variable]["viable"]
                data_dict[variable]["means"] = total ./ no_viable_simulations
            else
                data_dict[variable]["means"] = total ./ no_simulations
            end
        elseif data_dict[variable]["dims"] == 3
            means = zeros(no_reactions, size(total, 2))
            # Calculate means (first checking what should be divided by)
            if data_dict[variable]["divide_by_R"]
                for i in 1:no_reactions
                    means[i, :] = total[i, :] ./ no_simulations_with_R[i, :]
                end
            elseif data_dict[variable]["viable"]
                for i in 1:no_reactions
                    means[i, :] = total[i, :] ./ no_viable_simulations
                end
            else
                for i in 1:no_reactions
                    means[i, :] = total[i, :] ./ no_simulations
                end
            end
            # Add means to data dictionary
            data_dict[variable]["means"] = means
        end
    end
    return (data_dict)
end

# Written own function for standard deviation as it was tricky to get the standard function
# to handle no species simulations properly
function find_stand_dev(values::Vector{Float64}, mean::Float64, no_points::Float64)
    return sqrt(sum((values .- mean) .^ 2) / (no_points - 1))
end

# This function calculates the standard deviations for all variables along the trajectory
function calculate_trajectory_standard_devs!(data_dict::Dict, times::Vector{Float64},
                                             final_time_points::Vector{Float64},
                                             no_simulations::Vector{Float64},
                                             no_viable_simulations::Vector{Float64},
                                             no_simulations_with_R::Matrix{Float64},
                                             no_reactions::Int64)
    # Preallocate containers for the standard deviations
    for variable in keys(data_dict)
        data_dict[variable]["sds"] = zeros(size(data_dict[variable]["means"]))
    end
    # Loop over times
    for i in eachindex(times)
        # Find indices of still progressing trajectories
        inds = (final_time_points .>= times[i])
        # Find indices of still progressing trajectories with one or more viable strains
        vinds = (final_time_points .>= times[i]) .&
                (data_dict["viable_species"]["combined_data"][:, i] .> 0.0)
        # Find indices for where viable species in each reaction class exist
        rinds = Vector{Vector}(undef, no_reactions)
        for j in 1:no_reactions
            rinds[j] = (final_time_points .>= times[i]) .&
                       (data_dict["viable_species_per_reac_class"]["combined_data"][:, j,
                                                                                    i] .>
                        0.0)
        end
        # Loop over all variables
        for variable in keys(data_dict)
            # Calculate standard deviations (procedure changes based on what needs to be
            # averaged over)
            if data_dict[variable]["dims"] == 2 && data_dict[variable]["viable"]
                # These should be calculated just for viable strains
                if no_viable_simulations[i] > 1
                    sd_value = find_stand_dev(data_dict[variable]["combined_data"][vinds,
                                                                                   i],
                                              data_dict[variable]["means"][i],
                                              no_viable_simulations[i])
                    data_dict[variable]["sds"][i] = sd_value
                else
                    data_dict[variable]["sds"][i] = NaN
                end
            elseif data_dict[variable]["dims"] == 2
                sd_value = find_stand_dev(data_dict[variable]["combined_data"][inds, i],
                                          data_dict[variable]["means"][i],
                                          no_simulations[i])
                data_dict[variable]["sds"][i] = sd_value
            elseif data_dict[variable]["dims"] == 3 && data_dict[variable]["divide_by_R"]
                # Calculate standard deviations for reactions
                for j in 1:no_reactions
                    # Use only these in the reaction calculation
                    if no_simulations_with_R[j, i] > 1
                        sd_value = find_stand_dev(data_dict[variable]["combined_data"][rinds[j],
                                                                                       j,
                                                                                       i],
                                                  data_dict[variable]["means"][j, i],
                                                  no_simulations_with_R[j, i])
                        data_dict[variable]["sds"][j, i] = sd_value
                    else
                        data_dict[variable]["sds"][j, i] = NaN
                    end
                end
            elseif data_dict[variable]["dims"] == 3 && data_dict[variable]["viable"]
                # These should be calculated just for viable strains
                if no_viable_simulations[i] > 1
                    for j in 1:no_reactions
                        sd_value = find_stand_dev(data_dict[variable]["combined_data"][vinds,
                                                                                       j,
                                                                                       i],
                                                  data_dict[variable]["means"][j, i],
                                                  no_viable_simulations[i])
                        data_dict[variable]["sds"][j, i] = sd_value
                    end
                else
                    data_dict[variable]["sds"][:, i] .= NaN
                end
            elseif data_dict[variable]["dims"] == 3
                for j in 1:no_reactions
                    sd_value = find_stand_dev(data_dict[variable]["combined_data"][inds, j,
                                                                                   i],
                                              data_dict[variable]["means"][j, i],
                                              no_simulations[i])
                    data_dict[variable]["sds"][j, i] = sd_value
                end
            end
        end
    end
    return (data_dict)
end

# Function to read in variables with time and calculate stats over time
function calculate_trajectory_stats()
    # Check that sufficient arguments have been provided
    if length(ARGS) < 3
        error("insufficient inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    repeats = 0
    ims = 0
    sim_type = 0
    # Check that all arguments can be converted to integers
    try
        repeats = parse(Int64, ARGS[1])
        ims = parse(Int64, ARGS[2])
        sim_type = parse(Int64, ARGS[3])
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
    no_steps = 2500
    # Container to store final times
    final_time_points = zeros(repeats)
    # Define here so that it is available outside the for loop
    no_reactions = 0
    # Loop over number of repeats
    for i in 1:repeats
        # Load in relevant output file
        vfile = "Output/$(tk)$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/AvRun$(i)Data$(ims)Ims.jld"
        if ~isfile(vfile)
            error("$(ims) immigrations run $(i) is missing a variables file")
        end
        # Just want to save final times for now
        final_time_points[i] = load(vfile, "final_time_point")
        # Save number of reactions from the first run
        if i == 1
            viable_species_per_reac_class = load(vfile, "viable_species_per_reac_class")
            no_reactions = size(viable_species_per_reac_class, 1)
        end
    end
    # Use maximum final time to set time value
    times = collect(range(0.0, maximum(final_time_points), length = no_steps))
    # Preallocate relevant containers
    no_simulations = zeros(length(times))
    no_viable_simulations = zeros(length(times))
    no_simulations_with_R = zeros(no_reactions, length(times))
    # Define the other variables of interest + their required dimensionality
    # + whether to divide by only simulations with viable species
    # + whether to divide by number of species in reaction number class
    variables_of_interest = [
        ("surviving_species", 2, false, false),
        ("viable_species", 2, false, false),
        ("total_population", 2, false, false),
        ("shannon_diversity", 2, false, false),
        ("no_substrates", 2, false, false),
        ("species_per_reac_class", 3, false, false),
        ("viable_species_per_reac_class", 3, true, false),
        ("average_no_reac_steps", 2, true, false),
        ("average_η", 2, true, false),
        ("average_ω", 2, true, false),
        ("total_biomass_of_viable_species", 2, true, false),
        ("average_ΔG", 2, true, false),
        ("average_η_per_reac_class", 3, true, true),
        ("average_KS_per_reac_class", 3, true, true),
    ]
    # Convert this list into a dictionary of preallocated arrays
    data_dict = make_data_dictonary_from_list(variables_of_interest, repeats,
                                              no_reactions,
                                              length(times))
    # Final ϕR values are a special case
    all_final_ϕRs = Float64[]
    # Loop over number of trajectories (to minimise the number of reads in)
    for i in 1:repeats
        # Load in relevant output file
        vfile = "Output/$(tk)$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/AvRun$(i)Data$(ims)Ims.jld"
        if ~isfile(vfile)
            error("$(ims) immigrations run $(rN) is missing a variables file")
        end
        # First find and save final ϕR value for run
        final_ϕR = load(vfile, "final_ϕR")
        all_final_ϕRs = cat(all_final_ϕRs, final_ϕR, dims = 1)
        # Then load the time interval data
        T = load(vfile, "T")
        # Load all other variables into the dictionary
        data_dict = load_trajectory_vars_to_dict!(vfile, data_dict)
        # Bool to indicate end of the run
        run_end = false
        cnt = 0
        # Loop until end of run reached
        while run_end == false
            # increment counter
            cnt += 1
            # Still a trajectory here so increment that count by one
            no_simulations[cnt] += 1
            # Find index of first point greater than or equal to time
            Tind = findfirst(x -> x >= times[cnt], T)
            # Only increment the viable counter if there are viable strains at time point
            if data_dict["viable_species"]["run_data"][Tind] != 0
                no_viable_simulations[cnt] += 1
            end
            # Loop over reactions to find number of viable reactions
            for j in 1:no_reactions
                # Check if counter should be incremented
                if data_dict["viable_species_per_reac_class"]["run_data"][j, Tind] .> 0.0
                    no_simulations_with_R[j, cnt] += 1
                end
            end
            # Use function to add data into combined data container
            data_dict = add_to_combined_data!(data_dict, T, times, Tind, cnt, i)
            # Finally check if next time point is higher than final time for this trajectory
            if cnt >= length(times) || times[cnt + 1] > final_time_points[i]
                run_end = true
            end
        end
        println("Analysed trajectory $(i)")
        flush(stdout)
    end
    # Use function to calculate means
    data_dict = calculate_trajectory_means!(data_dict, no_reactions, no_simulations,
                                            no_viable_simulations, no_simulations_with_R)
    println("Means found")
    # Then use function to calculate standard deviations
    data_dict = calculate_trajectory_standard_devs!(data_dict, times, final_time_points,
                                                    no_simulations, no_viable_simulations,
                                                    no_simulations_with_R, no_reactions)
    println("Standard deviations found")
    # Now want to save means and standard deviations
    jldopen("Output/$(tk)$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/RunStats$(ims)Ims.jld",
            "w") do file
        # Save times
        write(file, "times", times)
        # Save number of continuing trajectories
        write(file, "no_simulations", no_simulations)
        write(file, "no_viable_simulations", no_viable_simulations)
        write(file, "no_simulations_with_R", no_simulations_with_R)
        # Save means and standard deviations
        for variable in keys(data_dict)
            write(file, "mean_$(variable)", data_dict[variable]["means"])
            write(file, "sd_$(variable)", data_dict[variable]["sds"])
        end
        # Finally write all of the final ϕR values out
        write(file, "all_final_ϕRs", all_final_ϕRs)
    end
    println("All data saved")
    return (nothing)
end

# Function to read in snapshot data and calculate stats over time
function snpstats()
    # Check that sufficient arguments have been provided
    if length(ARGS) < 3
        error("insufficient inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    repeats = 0
    ims = 0
    sim_type = 0
    # Check that all arguments can be converted to integers
    try
        repeats = parse(Int64, ARGS[1])
        ims = parse(Int64, ARGS[2])
        sim_type = parse(Int64, ARGS[3])
    catch e
        error("need to provide 3 integers")
    end
    println("Compiled")
    flush(stdout)
    # Load in hardcoded simulation parameters
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Load in 1st output file
    sfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/SnapData$(ims)Ims.jld"
    if ~isfile(sfile)
        error("$(ims) immigrations run 1 is missing a snapshot data file")
    end
    # Load in snapshot times
    times = load(sfile, "times")
    ns = load(sfile, "ns")
    gs = load(sfile, "gs")
    stb = load(sfile, "stb")
    inc = load(sfile, "inc")
    dec = load(sfile, "dec")
    st_r = load(sfile, "st_r")
    # Construct totals here
    tot_stb = dropdims(sum(stb, dims = 2), dims = 2)
    tot_inc = dropdims(sum(inc, dims = 2), dims = 2)
    tot_dec = dropdims(sum(dec, dims = 2), dims = 2)
    # Now calculate means
    mn_stb = tot_stb ./ st_r
    mn_inc = tot_inc ./ st_r
    mn_dec = tot_dec ./ st_r
    println("Means found")
    # Preallocate containers for the standard deviations
    sd_stb = zeros(size(mn_stb))
    sd_inc = zeros(size(mn_inc))
    sd_dec = zeros(size(mn_dec))
    # Loop over times
    for i in 1:(length(times) - 1)
        # Find indices of still progressing trajectories
        inds = (ns[i, :] .!== 0.0)
        # Calculate standard deviations
        sd_stb[i] = sqrt(sum((stb[i, inds] .- mn_stb[i]) .^ 2) / (st_r[i] - 1))
        sd_inc[i] = sqrt(sum((inc[i, inds] .- mn_inc[i]) .^ 2) / (st_r[i] - 1))
        sd_dec[i] = sqrt(sum((dec[i, inds] .- mn_dec[i]) .^ 2) / (st_r[i] - 1))
    end
    # Find total numbers of new strains and total number that grows
    tot_ns = dropdims(sum(ns, dims = 2), dims = 2)
    tot_gs = dropdims(sum(gs, dims = 2), dims = 2)
    # Use to find growth probability
    gp = zeros(size(tot_gs))
    for i in eachindex(gp)
        # Check that value isn't zero
        if tot_ns[i] != 0
            gp[i] = tot_gs[i] / tot_ns[i]
        else
            gp[i] = 0.0
        end
    end
    # Now just save the relevant data
    jldopen("Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/SnapDataStats$(ims)Ims.jld",
            "w") do file
        # Save times of snapshots
        write(file, "times", times)
        # Save growth probabilities
        write(file, "gp", gp)
        # Save means
        write(file, "mn_stb", mn_stb)
        write(file, "mn_inc", mn_inc)
        write(file, "mn_dec", mn_dec)
        # Save standard deviations
        write(file, "sd_stb", sd_stb)
        write(file, "sd_inc", sd_inc)
        write(file, "sd_dec", sd_dec)
        # Save number of simulations
        write(file, "st_r", st_r)
    end
    return (nothing)
end

@time calculate_trajectory_stats()
