# Script to read in organise and plot the ATP data
using CSV
using Plots
import PyPlot

function atp_read()
    # Read in csv file
    dataf = CSV.read("Data/dataset_22C_3d_atp_2.csv")
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
            rps = unique(dataf.replicate)
            nr = length(rps)
            # Then loop over the replicates
            for k = 1:nr
                # Find data for particular index
                I2 = I .& (dataf.replicate .== rps[k])
                # Then plot the data
                if i == 3
                    scatter!(dataf.minute[I2]/60.0,dataf.trait_value[I2]*1e-9,label="")
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
