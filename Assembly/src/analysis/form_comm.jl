# A script to form stable communities using the full model that then saves them for later use
using Assembly
using JLD
using SymPy

# Function to assemble specfic communities
function assemble()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 3
        error("need to specify number of reactions, number of strains and number of repeats")
    end
    # Preallocate the variables I want to extract from the input
    R = 0
    N = 0
    rps = 0
    # Check that all arguments can be converted to integers
    try
        R = parse(Int64,ARGS[1])
        N = parse(Int64,ARGS[2])
        rps = parse(Int64,ARGS[3])
    catch e
           error("all three inputs must be integer")
    end
    # Check that simulation type is valid
    if R < 1
        error("invalid number of reactions each strain must have more than 1 reaction")
    end
    # Check that number of strains is greater than 0
    if N < 1
        error("number of strains should be greater than zero")
    end
    # Check that number of strains is greater than 0
    if rps < 1
        error("need to do at least 1 simulation")
    end
    # Now start actual script
    println("Compiled and input read in!")
    flush(stdout)
    # Assume that half saturation occurs at a quarter κ/δ
    KS = (1/4)*5.5e-3
    # From wikipedia an average enzyme has a k to KS ratio of 10^5 M^-1 s^-1
    # This would give us a k of 137.5, sensible to assume an above average rate
    # Though should be reduced by the fact we include uptake as well as metabolism
    # Choosing k = 500 means we match the maximum glucose uptake rate seen in Natarajan et al (2000)
    # of 3*10^7 molecules per second.
    # The above is a sensible argument but 1.0 gives a more reasonable ATP concentration.
    kc = 10.0
    # The reversibility factor remains the same as previously
    kr = 10.0
    # Assume microbes have 3 reactions each
    # ALSO NEED TO WORK ON DIFFERENT CHOICES OF η
    mR = convert(Float64,R)
    sdR = 0.0
    # Case of 8 metabolites
    M = 8
    # Use formula to find how many reactions this implies
    O = 2*M - 3
    # Set time long enough for dynamics to equilbrate
    Tmax = 1e8
    # Fairly arbitary inital conditions
    pop = ones(N)
    conc = zeros(M)
    as = 1e5*ones(N)
    ϕs = 0.1*ones(N)
    # Now loop over the number of repeats
    for i = 1:rps
        # Print that the new run has been started
        println("Run $i started!")
        flush(stdout)
        # Make parameter set
        ps = initialise(N,M,O,mR,sdR,kc,KS,kr)
        # Before running the parameter sets should be saved so that if they crash
        # they can be rerun and hopefully track down where they went wrong
        # ONCE I'VE (HOPEFULLY) SOLVED THE PROBLEM DELETE THIS SECTION
        jldopen("Paras/ParasType$(R)Run$(i).jld","w") do file
            write(file,"ps",ps)
        end
        # Find starting time
        ti = time()
        # Then run the simulation
        C, T = full_simulate(ps,Tmax,pop,conc,as,ϕs)
        # And then print time elapsed
        tf = time()
        println("Time elapsed on run $i: $(tf-ti) s")
        flush(stdout)
        # Establish which microbes are extinct
        ext = (C[end,1:N] .== 0.0)
        # Preallocate vector to store extinct microbes
        ded = Array{MicrobeP,1}(undef,sum(ext))
        # Loop over and store microbes in the vector
        k = 0
        for j = 1:length(ext)
            if ext[j] == 1
                k += 1
                ded[k] = ps.mics[j]
            end
        end
        # Remove extinct strains from parameter set
        ps = extinction(ps,ext)
        # Preallocate final concentrations (etc) for output
        out = Array{Float64,1}(undef,3*ps.N+M)
        # Store final metabolite concentrations
        out[ps.N+1:ps.N+M] = C[end,N+1:N+M]
        # Now sub in data for not extinct microbes
        k = 0
        for j = 1:length(ext)
            if ext[j] != 1
                k += 1
                # Population
                out[k] = C[end,j]
                # Energy
                out[M+ps.N+k] = C[end,M+N+j]
                # Fraction
                out[M+2*ps.N+k] = C[end,M+2*N+j]
            end
        end
        # Save extinct strains
        jldopen("Output/ExtinctType$(R)Run$(i).jld","w") do file
            write(file,"ded",ded)
        end
        # the reduced parameter sets
        jldopen("Paras/ParasType$(R)Run$(i).jld","w") do file
            write(file,"ps",ps)
        end
        # and the full output
        jldopen("Output/OutputType$(R)Run$(i).jld","w") do file
            # Save final output
            write(file,"out",out)
            # # Save time data and dynamics data
            write(file,"T",T)
            write(file,"C",C[1:end,1:end])
        end
        # Print to show that run has been successfully completed
        println("Run $i completed and saved!")
        flush(stdout)
    end
    return(nothing)
end

# Function to find the stability times of the various simulations
function stability()
    println("Starting stability checking script")
    flush(stdout)
    # Check that sufficent arguments have been provided
    if length(ARGS) < 2
        error("need to specify community and number of repeats")
    end
    # Preallocate the variables I want to extract from the input
    R = 0
    rps = 0
    # Check that all arguments can be converted to integers
    try
        R = parse(Int64,ARGS[1])
        rps = parse(Int64,ARGS[2])
    catch e
           error("both inputs must be integer")
    end
    # Check that simulation type is valid
    if R < 1
        error("each strain must have more than 1 reaction")
    end
    # Check that number of simulations is greater than 0
    if rps < 1
        error("need to do at least 1 simulation")
    end
    # Preallocate stability times
    sT = zeros(rps)
    # Now loop over repeats
    for i = 1:rps
        # First check that files exists
        pfile = "Paras/ParasType$(R)Run$(i).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        ofile = "Output/OutputType$(R)Run$(i).jld"
        if ~isfile(ofile)
            error("run $(i) is missing an output file")
        end
        efile = "Output/ExtinctType$(R)Run$(i).jld"
        if ~isfile(ofile)
            error("run $(i) is missing an extinct file")
        end
        # Then load in data
        ps = load(pfile,"ps")
        out = load(ofile,"out")
        # Now preallocate vector of forces
        F = Array{Sym,1}(undef,3*ps.N+ps.M)
        # Find forces using function
        F = Force(ps,F)
        # Use final simulation results to find local forces
        f = nForce(F,out,ps)
        # Check maximum force on population
        stabN = (maximum(abs.(f[1:ps.N])) <= 2.0)
        # And the maximum force on the energies
        staba = (maximum(abs.(f[(ps.N+ps.M+1):(2*ps.N+ps.M)])) <= 1.0e-5)
        # Check if the populations and the energies are stable
        if staba == true && stabN == true
            # Find corresponding data
            T = load(ofile,"T")
            C = load(ofile,"C")
            # Setup while loop
            j = length(T)
            stab = false
            # Find orginal number of species
            N = round(Int64,(size(C)[2]-ps.M)/3)
            # Find index where final population matches saved out put
            inds = indexin(out[1:ps.N],C[end,1:N])
            println("Number of time points $(length(T))")
            flush(stdout)
            # Now want to track down the first time point that the system stabilises at
            # THIS NEEDS TO CHANGE TO RUN BACKWARDS FOR THE SAKE OF SPEED
            while stab == false
                j -= 1
                # Select the part of the data that we care about
                out = [C[j,inds];C[j,N+1:N+ps.M];C[j,inds.+N.+ps.M];C[j,inds.+2*N.+ps.M]]
                # If the above is wrong expect a "Wrong length" error here
                f = nForce(F,out,ps)
                # Check maximum force on population
                stabN = (maximum(abs.(f[1:ps.N])) <= 2.0)
                # And the maximum force on the energies
                staba = (maximum(abs.(f[(ps.N+ps.M+1):(2*ps.N+ps.M)])) <= 1.0e-5)
                # Looking for first unstable point
                if staba == false || stabN == false
                    stab = true
                    # Add previous time step to the list of stable times
                    sT[i] = T[j+1]
                    println("Time point for stability $(j+1)")
                    flush(stdout)
                end
            end
        else
            # For systems that arn't stable leave blank
            sT[i] = NaN
        end
    end
    # Calculate and output the number of stable systems
    Ns = length(sT)-sum(isnan.(sT))
    println("$(Ns) systems out of $(rps) are stable")
    # Output the data
    jldopen("Output/StabTimesType$(R).jld","w") do file
        write(file,"sT",sT)
    end
    return(nothing)
end

if length(ARGS) == 3
    @time assemble()
elseif length(ARGS) == 2
    @time stability()
end
