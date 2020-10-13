# Script to read in organise and plot the ATP data
using CSV
using DataFrames
using Plots
using Statistics
using StatsBase
import PyPlot

# Writing a function to perform a hampel filter on a data set
# IF I USE THIS IN A PAPER I NEED TO CHECK IT THROUGH CAREFULLY
function hampel_filter(x::Array{Float64,1},k::Int64)
    # Preallocate vector of points to keep
    kp = fill(true,length(x))
    # First and last point are preserved
    for i = 2:(length(x)-1)
        # Find start point of the sample
        if i <= k
            st = 1
        else
            st = i - k
        end
        # Find end point of the sample
        if i + k >= length(x)
            ed = length(x)
        else
            ed = i + k
        end
        # Now make sample
        sx = x[st:ed]
        # Find median of the sample
        mx = median(sx)
        # Find median absolute deviation of the sample about the median
        mdx = mad(x;center=mx,normalize=false)
        # Convert into standard deviation
        sdx = 1.4826*mdx
        # Check if point lies within 3 standard deviations of sample median
        if x[i] > mx + 3*sdx || x[i] < mx - 3*sdx
            kp[i] = false
        end
    end
    return(kp)
end

# Need to work out how to remove the correct points
function atp_read()
    # Read in csv file
    dataf = DataFrame!(CSV.File("Data/dataset_22C_3d_atp_2.csv"))
    # Count the number of unique ID's
    ids = unique(dataf.ID)
    ni = length(ids)
    # count the number of trait names
    tns = unique(dataf.trait_name)
    nt = length(tns)
    # Need to make a data structure to contain all the data
    data = fill(Float64[],ni,nt,4,2)
    # Structure to store genus names
    genus = Array{String,1}(undef,ni)
    # Store genus names
    for i = 1:ni
        # Find all data of a certain type for a particular ID
        It = (dataf.ID .== ids[i])
        # Find genus name to use as a title
        genus[i] = dataf.bacterial_genus[It][1]
    end
    # Make a bool of valid points
    val = fill(true,length(dataf.ID))
    # Locate anamolous points
    an1 = (dataf.bacterial_genus .== "Flavobacterium") .& (dataf.replicate .== 2) .&
    (dataf.minute .>= 3012) .& (dataf.trait_name .== tns[3])
    an2 = (dataf.bacterial_genus .== "Pseudomonas(2)") .& (dataf.replicate .== 1) .&
    (dataf.minute .>= 3012) .& (dataf.trait_name .== tns[3])
    an3 = (dataf.bacterial_genus .== "Serratia") .& (dataf.replicate .== 1) .&
    (dataf.trait_name .== tns[3])
    # Set these points as invalid
    val[an1] .= false
    val[an2] .= false
    val[an3] .= false
    # Loop over all three conditions
    for m = 1:3
        for i = 1:ni
            # Find all cell count data with species ID
            I = (dataf.ID .== ids[i]) .& (dataf.trait_name .== tns[m])
            # Count number of replicates
            rps = unique(dataf.replicate[I])
            nr = length(rps)
            # Find unique times and order them
            ts = sort(unique(dataf.minute[I]))
            # Only remove points if they are ATP data
            if m == 3
                # Now loop over times
                for j = 1:length(ts)
                    # Group samples by time
                    IT = I .& (dataf.minute .== ts[j])
                    # Then perform hampel filter on sample
                    kp = hampel_filter(dataf.trait_value[IT],3)
                    # Loop over this output
                    for k = 1:length(kp)
                        # Overwrite the replicate number with 0
                        if kp[k] == 0
                            # Find indices of true values in IT
                            inds = findall(x->x==true,IT)
                            # The k'th index should then be marked false
                            val[inds[k]] = false
                        end
                    end
                end
            end
            # Loop over replicates
            for j = 1:nr
                # Only save valid cases with particular replicate number
                I2 = I .& (dataf.replicate .== rps[j]) .& val
                # Save times
                data[i,m,j,1] = dataf.minute[I2]
                # and save corresponding values
                if m == 3
                    data[i,m,j,2] = 1e-9.*dataf.trait_value[I2]
                else
                    data[i,m,j,2] = dataf.trait_value[I2]
                end
            end
        end
    end
    # Call PyPlot
    pyplot()
    # Set a color-blind friendly palette
    theme(:wong2,dpi=150)
    # Plot all graphs
    for k = 1:3
        for i = 1:ni
            # Add appropriate times and labels
            plot(xlabel="Time (hours)",ylabel=tns[3],title=genus[i])
            # Loop over repeats to plot replicates
            for j = 1:4
                scatter!(data[i,k,j,1]/60.0,data[i,k,j,2],label="")
            end
            # Then save plot
            if k == 1
                savefig("Output/ATPDataPlots/Cells/$(genus[i]).png")
            elseif k == 2
                savefig("Output/ATPDataPlots/Biomass/$(genus[i]).png")
            else
                savefig("Output/ATPDataPlots/ATP/$(genus[i]).png")
            end
        end
    end
    return(nothing)
end

@time atp_read()
