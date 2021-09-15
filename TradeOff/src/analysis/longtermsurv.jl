# Script to find and save long term survival data
using TradeOff
using JLD

# Something
function long_term_surv()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 3
        error("insufficent inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    rps = 0
    ims = 0
    sim_type = 0
    # Check that all arguments can be converted to integers
    try
        rps = parse(Int64,ARGS[1])
        ims = parse(Int64,ARGS[2])
        sim_type = parse(Int64,ARGS[3])
    catch e
        error("need to provide 3 integers")
    end
    # Load in hardcoded simulation parameters
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    println("Compiled")
    flush(stdout)
    # Preallocate survival times
    sTs = Float64[]
    # Loop over number of repeats
    for i = 1:rps
        # Load in relevant output file
        ofile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Run$(i)Data$(ims)Ims.jld"
        if ~isfile(ofile)
            error("$(ims) immigrations run $(rN) is missing an output file")
        end
        # Load in microbe data
        micd = load(ofile,"micd")
        # Find indices of species that survive to the end
        inds = findall(isnan,micd.↦:ExT)
        # Add the immigration times of these species to the vector
        sTs = cat(sTs,(micd.↦:ImT)[inds],dims=1)
    end
    # Once all this has been calculated save ``survival'' times as a new datafiles
    # Now want to save means and standard deviations
    jldopen("Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/SurvTimes$(ims)Ims.jld","w") do file
        # Save times
        write(file,"sTs",sTs)
    end
    return(nothing)
end

@time long_term_surv()
