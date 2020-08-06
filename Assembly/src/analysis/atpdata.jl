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

function atp_read()
    # Read in csv file
    # dataf = CSV.read("Data/dataset_22C_3d_atp_2.csv")
    dataf = DataFrame!(CSV.File("Data/dataset_22C_3d_atp_2.csv"))
    # Count the number of unique ID's
    ids = unique(dataf.ID)
    ni = length(ids)
    # count the number of trait names
    tns = unique(dataf.trait_name)
    nt = length(tns)
    # Call PyPlot
    pyplot(dpi=150)
    # Then plot all the graphs at once
    for i = 1:nt
        for j = 1:ni
            # Find all data of a certain type for a particular ID
            I = (dataf.ID .== ids[j]) .& (dataf.trait_name .== tns[i])
            # Find genus name to use as a title
            gen = dataf.bacterial_genus[I][1]
            # Find trait name to use as a y label
            yl = dataf.trait_name[I][1]
            # Add the appropriate titles
            plot(xlabel="Time (hours)",ylabel=yl,title=gen)
            # Count number of replicates
            rps = unique(dataf.replicate[I])
            nr = length(rps)
            # Find unique times and order them
            ts = sort(unique(dataf.minute[I]))
            # Now loop over times
            for k = 1:length(ts)
                # Group samples by time
                IT = I .& (dataf.minute .== ts[k])
                # Then perform hampel filter on sample
                kp = hampel_filter(dataf.trait_value[IT],3)
                # Loop over this output
                for m = 1:length(kp)
                    # Overwrite the replicate number with 0
                    if kp[m] == 0
                        (dataf.replicate[IT])[m] = 0
                    end
                end
            end
            # IT DOESN'T FULLY WORK, I THINK THIS IS BECAUSE I NEED TO RUN IT ON
            # ONLY THE CELL COUNT AND THEN DELETE CORRESPONDING ATP VALUES
            # Now loop over the replicates
            for k = 1:nr
                # Find data for particular index
                I2 = I .& (dataf.replicate .== rps[k])
                # Then plot the data
                if i == 3
                    scatter!((dataf.minute[I2])/60.0,(dataf.trait_value[I2])*1e-9,label="")
                elseif i == 1
                    scatter!((dataf.minute[I2])/60.0,(dataf.trait_value[I2]),label="")
                else
                    scatter!(dataf.minute[I2]/60.0,dataf.trait_value[I2],label="")
                end
            end
            # Finally save the plot
            if i == 1
                savefig("Output/ATPDataPlots/$(gen)Cells.png")
            elseif i == 2
                savefig("Output/ATPDataPlots/$(gen)Biomass.png")
            else
                savefig("Output/ATPDataPlots/$(gen)ATP.png")
            end
        end
    end
    return(nothing)
end

@time atp_read()
