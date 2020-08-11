# A script to form stable communities using the full model that then saves them for later use
using Assembly
using Plots
using LaTeXStrings
using JLD
import PyPlot

# function to test that the new stuff I'm writing actually works
# Keeping this for the moment as the plotting stuff might be useful
function test()
    println("Compiled!")
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
    # Assume microbes have 2 reactions each
    mR = 2.0
    sdR = 0.0
    # Case of 8 metabolites and 1 strain
    N = 3
    M = 8
    # Use formula to find how many reactions this implies
    O = 2*M - 3
    good = false
    # Make parameter set
    ps = 0
    while good == false
        ps = initialise(N,M,O,mR,sdR,kc,KS,kr)
        # This step is to ensure that the first metabolite can be broken down
        for i = 1:N
            if any((ps.reacs[ps.mics[i].Reacs].↦:Rct) .== 1)
                good = true
            end
        end
    end
    # Set time long enough to see dynamics
    Tmax = 10000000.0
    # Fairly arbitary inital conditions
    pop = ones(N)
    conc = zeros(M)
    as = 1e5*ones(N)
    ϕs = 0.1*ones(N)
    C, T = full_simulate(ps,Tmax,pop,conc,as,ϕs)
    pyplot(dpi=200)
    # Setup population plot
    p1 = plot(xlabel="Time",ylabel="Population",yaxis=:log10)
    for i = 1:N
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0)
        plot!(p1,T[inds],C[inds,i],label="")
    end
    savefig(p1,"Output/PopvsTime.png")
    plot(T,C[:,N+1:N+M],xlabel="Time",label="",ylabel="Concentration")
    savefig("Output/MetabolitevsTime.png")
    plot(T,C[:,N+M+1:2*N+M],xlabel="Time",label="",ylabel="Cell energy conc")
    savefig("Output/EnergyvsTime.png")
    plot(T,C[:,2*N+M+1:end],xlabel="Time",label="",ylabel=L"\phi_R")
    savefig("Output/FractionvsTime.png")
    # Calculate growth rates
    λ = zeros(length(C[:,1]),N)
    for i = 1:length(C[:,1])
        for j = 1:N
            λ[i,j] = λs(C[i,N+M+j],C[i,2*N+M+j],ps.mics[j])
        end
    end
    # Then plot growth rates
    plot(T,λ,xlabel="Time",label="",ylabel=L"\lambda")
    savefig("Output/GrowthvsTime.png")
    # Calculate reaction inhibition
    # Find surviving microbes
    ms = []
    for i = 1:N
        if C[end,i] > 1e-8
            ms = cat(ms,i,dims=1)
        end
    end
    # Make array of bools to store results
    # To generalise this I would have to figure out how to make
    # an array with rows of differing lengths
    rs = fill(false,(length(ms),ps.mics[1].R))
    # Then for these microbes find viable reactions
    for i = 1:length(ms)
        # Find corresponding microbes
        mic = ps.mics[ms[i]]
        for j = 1:mic.R
            # find reaction
            r = ps.reacs[mic.Reacs[j]]
            if C[end,N+r.Rct] > 1e-8 && C[end,N+r.Prd] > 1e-8
                rs[i,j] = true
            end
        end
    end
    # Preallocate array of θs, this also will take some work to generalise
    θs = zeros(length(C[:,1]),size(rs,1)*size(rs,2))
    for i = 1:length(ms)
        for j = 1:ps.mics[1].R
            # Find appropriate reaction
            r = ps.reacs[ps.mics[ms[i]].Reacs[j]]
            for k = 1:length(C[:,1])
                θs[k,(i-1)*ps.mics[1].R+j] = θ_smooth(C[k,N+r.Rct],C[k,N+r.Prd],ps.T,ps.mics[ms[i]].η[j],r.ΔG0)
            end
        end
    end
    # Then plot inhibitions
    plot(T,θs,xlabel="Time",label="",ylabel=L"\theta")
    savefig("Output/InhibvsTime.png")
    return(nothing)
end

# Function to assemble specfic communities
function assemble()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 3
        error("need to specify type of community to assemble, number of species and number of repeats")
    end
    # Preallocate the variables I want to extract from the input
    st = 0
    N = 0
    rps = 0
    # Check that all arguments can be converted to integers
    try
        st = parse(Int64,ARGS[1])
        N = parse(Int64,ARGS[2])
        rps = parse(Int64,ARGS[3])
    catch e
           error("all three inputs must be integer")
    end
    # Check that simulation type is valid
    if st > 4 || st < 1
        error("invalid simulation type, specify 1 for low η, 2 for moderate η, 3 for high η, or 4 for mixed")
    end
    # Check that number of strains is greater than 0
    if N < 1
        error("number of strains should be greater than zero")
    end
    # Check that number of strains is greater than 0
    if rps < 1
        error("need to do at least 1 simulation")
    end
    # Now
    println("Compiled and input read in!")
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
    mR = 3.0
    sdR = 0.0
    # Case of 8 metabolites
    M = 8
    # Use formula to find how many reactions this implies
    O = 2*M - 3
    # Set time long enough for dynamics to equilbrate
    Tmax = 10000000.0
    # Fairly arbitary inital conditions
    pop = ones(N)
    conc = zeros(M)
    as = 1e5*ones(N)
    ϕs = 0.1*ones(N)
    # Now loop over the number of repeats
    for i = 1:rps
        # Make parameter set
        ps = initialise(N,M,O,mR,sdR,kc,KS,kr,st)
        # Then run the simulation
        C, T = full_simulate(ps,Tmax,pop,conc,as,ϕs)
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
        jldopen("Output/ExtinctType$(st)Run$(i).jld","w") do file
            write(file,"ded",ded)
        end
        # the reduced parameter sets
        jldopen("Paras/ParasType$(st)Run$(i).jld","w") do file
            write(file,"ps",ps)
        end
        # and the final output
        jldopen("Output/OutputType$(st)Run$(i).jld","w") do file
            write(file,"out",out)
        end
    end
    return(nothing)
end

# A function to read in saved data and interpret it
function interpret()
    println("Compiled!")
    # Do this by hand for the moment
    ps = load("Paras/ParasType1Run1.jld","ps")
    out = load("Output/OutputType1Run1.jld","out")
    println("Loaded data")
    # Use these to run a simulation
    Tmax = 10000.0
    pop = out[1:ps.N]
    conc = out[(ps.N+1):(ps.N+ps.M)]
    as = out[(ps.N+ps.M+1):(2*ps.N+ps.M)]
    ϕs = out[(2*ps.N+ps.M+1):end]
    C, T = full_simulate(ps,Tmax,pop,conc,as,ϕs)
    # Then plot the final output
    pyplot(dpi=200)
    # Setup population plot
    p1 = plot(xlabel="Time",ylabel="Population",yaxis=:log10)
    for i = 1:ps.N
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0)
        plot!(p1,T[inds],C[inds,i],label="")
    end
    savefig(p1,"Output/PopvsTime.png")
    plot(T,C[:,ps.N+1:ps.N+ps.M],xlabel="Time",label="",ylabel="Concentration")
    savefig("Output/MetabolitevsTime.png")
    plot(T,C[:,ps.N+ps.M+1:2*ps.N+ps.M],xlabel="Time",label="",ylabel="Cell energy conc")
    savefig("Output/EnergyvsTime.png")
    plot(T,C[:,2*ps.N+ps.M+1:end],xlabel="Time",label="",ylabel=L"\phi_R")
    savefig("Output/FractionvsTime.png")
    # Calculate growth rates
    λ = zeros(length(C[:,1]),ps.N)
    for i = 1:length(C[:,1])
        for j = 1:ps.N
            λ[i,j] = λs(C[i,ps.N+ps.M+j],C[i,2*ps.N+ps.M+j],ps.mics[j])
        end
    end
    # Then plot growth rates
    plot(T,λ,xlabel="Time",label="",ylabel=L"\lambda")
    savefig("Output/GrowthvsTime.png")
    return(nothing)
end

@time assemble()
