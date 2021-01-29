# A script to form stable communities using the full model that then saves them for later use
using Assembly
using JLD

# Function to assemble specfic communities
function assemble()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 5
        error("Insufficent inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    Rl = 0
    Ru = 0
    N = 0
    rps = 0
    syn = false
    # Check that all arguments can be converted to integers
    try
        Rl = parse(Int64,ARGS[1])
        Ru = parse(Int64,ARGS[2])
        N = parse(Int64,ARGS[3])
        rps = parse(Int64,ARGS[4])
        syn = parse(Bool,ARGS[5])
    catch e
           error("need to provide 4 integers and a bool")
    end
    # Starting run assumed to be 1
    Rs = 1
    # Check if run to start from has been provided
    if length(ARGS) > 5
        # Check this argument is an integer
        try
            Rs = parse(Int64,ARGS[6])
        catch e
            error("intial run number must be integer")
        end
    end
    # Check that simulation type is valid
    if Rl < 1
        error("lower bound on the number of reactions must be greater than 1")
    end
    # Check that simulation type is valid
    if Ru < Rl
        error("upper bound on the number of reactions can't be smaller than the lower")
    end
    # Check that number of strains is greater than 0
    if N < 1
        error("number of strains should be greater than zero")
    end
    # Check that number of strains is greater than 0
    if rps < 1
        error("need to do at least 1 simulation")
    end
    if Rs > rps
        error("starting run can't be higher than final run")
    end
    println("Reaction range = $(Rl)-$(Ru)")
    println("Syntrophy on = $(syn)")
    # Now start actual script
    println("Compiled and input read in!")
    flush(stdout)
    # Assume that half saturation occurs at a quarter κ/δ
    KS = (1/4)*5.5e-3
    # Arbitary number that seems to give decent survival
    kc = 10.0
    # The reversibility factor remains the same as previously
    kr = 10.0
    # Assume microbes have 3 reactions each
    # Case of 8 metabolites
    M = 25
    # Use formula to find how many reactions this implies
    O = 2*M - 3
    # Set time long enough for dynamics to equilbrate
    Tmax = 1e8
    # Initial ribosome fraction is taken from my ATP fits
    ϕR0 = 0.128
    # Fairly arbitary inital conditions
    pop = ones(N)
    conc = zeros(M)
    as = 1e5*ones(N)
    ϕs = ϕR0*ones(N)
    # Now loop over the number of repeats
    for i = Rs:rps
        # Print that the new run has been started
        println("Run $i started!")
        flush(stdout)
        # Make parameter set
        ps = initialise(N,M,O,Rl,Ru,kc,KS,kr,syn)
        # Before running the parameter sets should be saved so that if they crash
        # they can be rerun and hopefully track down where they went wrong
        jldopen("Paras/ParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(rps).jld","w") do file
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
        jldopen("Output/ExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(rps).jld","w") do file
            write(file,"ded",ded)
        end
        # the reduced parameter sets
        jldopen("Paras/ParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(rps).jld","w") do file
            write(file,"ps",ps)
        end
        # and the full output
        jldopen("Output/OutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(rps).jld","w") do file
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

@time assemble()
