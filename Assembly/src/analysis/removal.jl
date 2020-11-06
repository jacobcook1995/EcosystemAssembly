# Script to remove strains that don't survive to infinity from the data sets
using Assembly
using JLD
using SymPy

# NEEDS TO BE CONVERTED FOR NEW OUTPUT
# ALSO NEEDS TO SAVE DATA AT INFINITY!!!
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
            # Write out old data if stable
            # Save extinct strains
            jldopen("Data/Type$(R)/RedExtinctType$(R)Run$(i).jld","w") do file
                write(file,"ded",ded)
            end
            # the reduced parameter sets
            jldopen("Data/Type$(R)/RedParasType$(R)Run$(i).jld","w") do file
                write(file,"ps",ps)
            end
            # and the full output
            jldopen("Data/Type$(R)/RedOutputType$(R)Run$(i).jld","w") do file
                # Save final output
                write(file,"out",out)
                # # Save time data and dynamics data
                write(file,"T",T)
                write(file,"C",C[1:end,1:end])
            end
        else
            # Set a high final time
            Tmax = 100*maximum(T)
            # Store intial numbers of strains and metabolites
            N = ps.N
            M = ps.M
            # Extract initial conditions
            pop = out[1:ps.N]
            conc = out[(ps.N+1):(ps.N+ps.M)]
            as = out[(ps.N+ps.M+1):(2*ps.N+ps.M)]
            ϕs = out[(2*ps.N+ps.M+1):end]
            # Then run the simulation
            Cl, Tl = full_simulate(ps,Tmax,pop,conc,as,ϕs)
            # Establish which microbes are now extinct
            ext = (Cl[end,1:N] .== 0.0)
            # Preallocate vector to store extinct microbes
            ded2 = Array{MicrobeP,1}(undef,sum(ext))
            # Loop over and store microbes in the vector
            k = 0
            for j = 1:length(ext)
                if ext[j] == 1
                    k += 1
                    ded2[k] = ps.mics[j]
                end
            end
            # Remove extinct strains from parameter set
            ps = extinction(ps,ext)
            # Preallocate final concentrations (etc) for output
            nout = Array{Float64,1}(undef,3*ps.N+M)
            # Store final metabolite concentrations
            nout[ps.N+1:ps.N+M] = out[N+1:N+M]
            # Now sub in data for not extinct microbes
            k = 0
            for j = 1:length(ext)
                if ext[j] != 1
                    k += 1
                    # Population
                    nout[k] = out[j]
                    # Energy
                    nout[M+ps.N+k] = out[M+N+j]
                    # Fraction
                    nout[M+2*ps.N+k] = out[M+2*N+j]
                end
            end
            # Gather and output new reduceded data
            ded = cat(ded,ded2,dims=1)
            # Save extinct strains
            jldopen("Data/Type$(R)/RedExtinctType$(R)Run$(i).jld","w") do file
                write(file,"ded",ded)
            end
            # the reduced parameter sets
            jldopen("Data/Type$(R)/RedParasType$(R)Run$(i).jld","w") do file
                write(file,"ps",ps)
            end
            # and the full output
            jldopen("Data/Type$(R)/RedOutputType$(R)Run$(i).jld","w") do file
                # Save final output
                write(file,"out",nout)
                # # Save time data and dynamics data
                write(file,"T",T)
                write(file,"C",C[1:end,1:end])
            end
        end
    end
    return(nothing)
end

@time removal()
