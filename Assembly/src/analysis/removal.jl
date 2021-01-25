# Script to remove strains that don't survive to infinity from the data sets
using Assembly
using JLD
using SymPy

# function to read in data set and remove non-long term suvivors
function removal()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 4
        error("Insufficent inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    Rl = 0
    Ru = 0
    syn = false
    nR = 0
    # Check that all arguments can be converted to integers
    try
        Rl = parse(Int64,ARGS[1])
        Ru = parse(Int64,ARGS[2])
        syn = parse(Bool,ARGS[3])
        nR = parse(Int64,ARGS[4])
    catch e
           error("should provide 3 integers and a bool")
    end
    # Check that simulation type is valid
    if Rl < 1
        error("lower bound on reactions must be positive")
    end
    # Check that simulation type is valid
    if Ru < Rl
        error("lower bound can't be greater than upper bound")
    end
    # Check that number of simulations is greater than 0
    if nR < 1
        error("number of repeats cannot be less than 1")
    end
    println("Compiled!")
    flush(stdout)
    # Setup counter
    cnt = 0
    # Loop over repeats
    for i = 1:nR
        # Assume that output files don't already exist
        outp = false
        # Three output files to check the existence of
        outf1 = "Data/$(Rl)-$(Ru)$(syn)/RedExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        outf2 = "Data/$(Rl)-$(Ru)$(syn)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        outf3 = "Data/$(Rl)-$(Ru)$(syn)/RedOutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        # Check if all three exist
        if isfile(outf1) && isfile(outf2) && isfile(outf3)
            outp = true
        end
        # Print out every time parameter sets
        if i % 10 == 0
            println("Reached run $(i)")
            flush(stdout)
        end
        # Read in relevant files
        pfile = "Data/$(Rl)-$(Ru)$(syn)/ParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        ofile = "Data/$(Rl)-$(Ru)$(syn)/OutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        if ~isfile(ofile)
            error("run $(i) is missing an output file")
        end
        efile = "Data/$(Rl)-$(Ru)$(syn)/ExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
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
        # Skip further evaluation if output already exists
        if outp == true
            println("Simulation $(i) already has output")
            flush(stdout)
        elseif stab == true
            # Increment counter
            cnt += 1
            # Write out old data if stable
            # Save extinct strains
            jldopen("Data/$(Rl)-$(Ru)$(syn)/RedExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld","w") do file
                write(file,"ded",ded)
            end
            # the reduced parameter sets
            jldopen("Data/$(Rl)-$(Ru)$(syn)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld","w") do file
                write(file,"ps",ps)
            end
            # and the full output
            jldopen("Data/$(Rl)-$(Ru)$(syn)/RedOutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld","w") do file
                # Save final output
                write(file,"out",out)
                # This output is basically the output at infinity
                write(file,"inf_out",out)
                # Save time data and dynamics data
                write(file,"T",T)
                write(file,"C",C[1:end,1:end])
            end
        else
            println("Simulation $(i) unstable")
            flush(stdout)
            # Set a high final time
            Tmax = 10*maximum(T)
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
            # Remove any microbes below threshold
            for j = 1:ps.N
                if Cl[end,j] < 1e-5
                    Cl[end,j] = 0.0
                end
            end
            # Preallocate force vector
            F = Array{Sym,1}(undef,3*ps.N+ps.M)
            # Find forces using function
            F = Force(ps,F)
            # Use final condition to find local forces
            f = nForce(F,Cl[end,:],ps)
            # Check if ecosystem is now stable
            stab2 = true
            for j = 1:ps.N
                if Cl[end,j] > 0.0 && abs(f[j]/Cl[end,j]) > 1e-9
                    # Find magnitude of the change
                    mg = abs(Cl[end,j] - Cl[1,j])
                    # Check if the magnitude of change is greater than the minimum value
                    if mg > minimum(Cl[:,j])
                        println("Strain $(j) unstable")
                        println("Pop = $(Cl[end,j])")
                        flush(stdout)
                        stab2 = false
                    end
                end
            end
            # Set up counter
            c = 0
            # If not yet stable run simulation again
            while stab2 == false && c < 50
                # Increment counter
                c += 1
                println("Simulation $(i) being repeated time $(c)")
                flush(stdout)
                # Extract initial conditions
                pop = Cl[end,1:ps.N]
                conc = Cl[end,(ps.N+1):(ps.N+ps.M)]
                as = Cl[end,(ps.N+ps.M+1):(2*ps.N+ps.M)]
                ϕs = Cl[end,(2*ps.N+ps.M+1):end]
                # Then run the simulation
                Cl, Tl = full_simulate(ps,Tmax,pop,conc,as,ϕs)
                # Remove any microbes below threshold
                for j = 1:ps.N
                    if Cl[end,j] < 1e-5
                        Cl[end,j] = 0.0
                    end
                end
                # Preallocate force vector
                F = Array{Sym,1}(undef,3*ps.N+ps.M)
                # Find forces using function
                F = Force(ps,F)
                # Use final condition to find local forces
                f = nForce(F,Cl[end,:],ps)
                # Check if ecosystem is now stable
                stab2 = true
                for j = 1:ps.N
                    if Cl[end,j] > 0.0 && abs(f[j]/Cl[end,j]) > 1e-9
                        # Find magnitude of the change
                        mg = abs(Cl[end,j] - Cl[1,j])
                        # Check if the magnitude of change is greater than the minimum value
                        if mg > minimum(Cl[:,j])
                            println("Strain $(j) unstable")
                            println("Pop = $(Cl[end,j])")
                            flush(stdout)
                            stab2 = false
                        end
                    end
                end
            end
            if c == 50
                println("Run $(i) is possibly oscilliatory")
                flush(stdout)
            end
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
            # Preallocate output at infinity
            inf_out = Array{Float64,1}(undef,3*ps.N+ps.M)
            # Store final metabolite concentrations
            inf_out[ps.N+1:ps.N+M] = Cl[end,N+1:N+M]
            # Now sub in data for not extinct microbes
            k = 0
            for j = 1:length(ext)
                if ext[j] != 1
                    k += 1
                    # Population
                    inf_out[k] = Cl[end,j]
                    # Energy
                    inf_out[ps.M+ps.N+k] = Cl[end,M+N+j]
                    # Fraction
                    inf_out[ps.M+2*ps.N+k] = Cl[end,M+2*N+j]
                end
            end
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
            jldopen("Data/$(Rl)-$(Ru)$(syn)/RedExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld","w") do file
                write(file,"ded",ded)
            end
            # the reduced parameter sets
            jldopen("Data/$(Rl)-$(Ru)$(syn)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld","w") do file
                write(file,"ps",ps)
            end
            # and the full output
            jldopen("Data/$(Rl)-$(Ru)$(syn)/RedOutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld","w") do file
                # Save final output
                write(file,"out",nout)
                # Save the output at infinity here
                write(file,"inf_out",inf_out)
                # Save time data and dynamics data
                write(file,"T",T)
                write(file,"C",C[1:end,1:end])
            end
        end
    end
    println("$(cnt) out of $(nR) were already stable")
    flush(stdout)
    return(nothing)
end

@time removal()
