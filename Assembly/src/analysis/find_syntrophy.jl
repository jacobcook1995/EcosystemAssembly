# This script exists to investigate syntrophy in the simulated communities
using Assembly
using JLD
using Plots
import PyPlot

# A function to search for syntrophic pairs in the data
function find_syn()
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
    # Now move onto plotting
    pyplot()
    theme(:wong2,dpi=200)
    # Loop over repeats
    for i = 1:nR
        # Read in relevant files
        pfile = "Data/Type$(R)/RedParasType$(R)Run$(i).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        ofile = "Data/Type$(R)/RedOutputType$(R)Run$(i).jld"
        if ~isfile(ofile)
            error("run $(i) is missing an output file")
        end
        efile = "Data/Type$(R)/RedExtinctType$(R)Run$(i).jld"
        if ~isfile(efile)
            error("run $(i) is missing an extinct file")
        end
        # Basically just loading everything out as I'm not sure what I'll need
        ps = load(pfile,"ps")
        C = load(ofile,"C")
        T = load(ofile,"T")
        out = load(ofile,"out")
        ded = load(efile,"ded")
        # Find orginal number of strains
        N = ps.N + length(ded)
        # Vector to indicate obligate syntrophically consuming strains
        cons = fill(false,ps.N)
        # Loop over strains
        for j = 1:ps.N
            # Extract strain of interest
            m = ps.mics[j]
            # Loop over its reactions
            for k = 1:m.R
                # Find reaction
                r = ps.reacs[m.Reacs[k]]
                # Find θ values
                θr = θ(out[ps.N+r.Rct],out[ps.N+r.Prd],ps.T,m.η[k],r.ΔG0)
                # Check if theta is high enough for syntrophy to be a meaningful effect
                if θr > 5e-4
                    # Loop over all other strains
                    for l = 1:ps.N
                        # Check if any strain other than itself consume the microbe
                        if l != j && any((ps.reacs[ps.mics[l].Reacs].↦:Rct) .== r.Prd)
                            # Update consumer strains to match
                            cons[l] = true
                        end
                    end
                end
            end
            # Check if there are any consumer strains
            if any(cons .== true)
                # Set up plotting of old dynamics
                plot(yaxis=:log10)
                Tmax = 5e6
                # Extract initial conditions
                pop = out[1:ps.N]
                conc = out[(ps.N+1):(ps.N+ps.M)]
                as = out[(ps.N+ps.M+1):(2*ps.N+ps.M)]
                ϕs = out[(2*ps.N+ps.M+1):end]
                # Find indices of down stream
                println("Run $i")
                println(cons)
                # cons always comes out as [1,0,...,0]
                # SUSPIOUS
                inds = collect(1:ps.N)[cons]
                # Simulate with strain repressed
                Cn, Tn = full_simulate_syn(ps,Tmax,pop,conc,as,ϕs,inds)
                # Set offset time
                Toff = 1e6
                # and use to join old data to new data
                Ca = cat(out',out',Cn[:,:],dims=1)
                Ta = cat(0.0,Toff,Tn.+Toff,dims=1)
                # Setup plot for new dynamics
                plot(yaxis=:log10,xlabel="Time (seconds)",ylabel="Log population")
                for i = 1:ps.N
                    # Find and eliminate zeros so that they can be plotted on a log plot
                    inds = (Ca[:,i] .> 0)
                    plot!(Ta[inds],Ca[inds,i],label="")
                end
                vline!([Toff],color=:red,label="Extinction")
                savefig("Output/Newpops$(i).png")
                # NEED A WAY TO CHECK IF TRUELY STABLE
                # AND THEN NEED TO SEE IF θ VALUE CORRELATES WITH STABILITY
                # NEED TO CHECK THAT THE ABOVE RESULTS ACTUALLY MAKE SENSE
                # BIG TASK TOMORROW IS TO THINK ABOUT MULTI-REACTION CASE
            end
        end
    end
    return(nothing)
end

@time find_syn()
