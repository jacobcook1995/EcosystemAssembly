# Script to remove strains that don't survive to infinity from the data sets
using Assembly
using JLD
using SymPy
# REMOVE THESE WHEN TESTING IS DONE
using Plots
import PyPlot

# function to read in data set and remove non-long term suvivors
function removal()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 2
        error("need to specify community and number of repeats")
    end
    # Preallocate the variables I want to extract from the input
    R = 0
    nR = 0
    # Check that all arguments can be converted to integers
    try
        R = parse(Int64,ARGS[1])
        nR = parse(Int64,ARGS[2])
    catch e
           error("both inputs must be integer")
    end
    # Check that simulation type is valid
    if R < 1
        error("each strain must have more than 1 reaction")
    end
    # Check that number of simulations is greater than 0
    if nR < 1
        error("Number of repeats cannot be less than 1")
    end
    println("Compiled!")
    # Set up plotting, REMOVE WHEN DONE
    pyplot()
    # Set a color-blind friendly palette
    theme(:wong2,dpi=200)
    # Loop over repeats
    for i = 1:nR
        # Read in relevant files
        pfile = "Data/Type$(R)/ParasType$(R)Run$(i).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        ofile = "Data/Type$(R)/OutputType$(R)Run$(i).jld"
        if ~isfile(ofile)
            error("run $(i) is missing an output file")
        end
        efile = "Data/Type$(R)/ExtinctType$(R)Run$(i).jld"
        if ~isfile(efile)
            error("run $(i) is missing an extinct file")
        end
        # Basically just loading everything out as I'm not sure what I'll need
        ps = load(pfile,"ps")
        C = load(ofile,"C")
        T = load(ofile,"T")
        out = load(ofile,"out")
        ded = load(efile,"ded")
        # Preallocate force vector
        F = Array{Sym,1}(undef,3*ps.N+ps.M)
        # Find forces using function
        F = Force(ps,F)
        # Use initial conditions to find local forces
        f = nForce(F,out,ps)
        # Check if forces are stable
        stab = all(abs.(f[1:ps.N]./out[1:ps.N]) .< 1e-9)
        if stab == true
            
        else

        end
        # # Now section to write this new data out
        # # Save extinct strains
        # jldopen("Data/Type$(R)/RedExtinctType$(R)Run$(i).jld","w") do file
        #     write(file,"ded",ded)
        # end
        # # the reduced parameter sets
        # jldopen("Data/Type$(R)/RedParasType$(R)Run$(i).jld","w") do file
        #     write(file,"ps",ps)
        # end
    end
    return(nothing)
end

@time removal()
