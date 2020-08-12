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
    # Now start actual script
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

# function to calculate the dissipation for an assembled ecosystem
function dissipation(ps::FullParameters,out::Array{Float64,1})
    # check that parameter set is sensible given the output
    if length(out) != ps.M + 3*ps.N
        error("parameter set doesn't match out put")
    end
    # Set dissipation to zero
    dsp = 0
    # Loop over number of strains
    for i = 1:ps.N
        # Isolate this strain
        mic = ps.mics[i]
        # Loop over reactions of this strain
        for j = 1:mic.R
            # Find appropriate reaction
            r = ps.reacs[mic.Reacs[j]]
            # Find amount of energy that this reaction dissipates
            Fd = -(r.ΔG0 + Rgas*ps.T*log(out[ps.N+r.Prd]/out[ps.N+r.Rct]) + mic.η[j]*ΔGATP)
            # Find amount of enzyme E
            E = Eα(out[2*ps.N+ps.M+i],mic,j)
            # Then find the rate that this reaction proceeds at
            q = qs(out[ps.N+r.Rct],out[ps.N+r.Prd],E,j,mic,ps.T,r)
            if q != 0.0
                dsp += q*Fd
            end
        end
    end
    return(dsp)
end

# A function to read in saved data and interpret it
function interpret()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 2
        error("need to specify community and number of repeats")
    end
    # Preallocate the variables I want to extract from the input
    st = 0
    rps = 0
    # Check that all arguments can be converted to integers
    try
        st = parse(Int64,ARGS[1])
        rps = parse(Int64,ARGS[2])
    catch e
           error("both inputs must be integer")
    end
    # Check that simulation type is valid
    if st > 4 || st < 1
        error("invalid simulation type, specify 1 for low η, 2 for moderate η, 3 for high η, or 4 for mixed")
    end
    # Check that number of simulations is greater than 0
    if rps < 1
        error("need to do at least 1 simulation")
    end
    println("Compiled!")
    # Preallocate memory to store the number of survivors in
    srv = zeros(rps)
    # Preallocate name
    nme = ""
    # Then find name for the community
    if st == 1
        nme = "Low η"
    elseif st == 2
        nme = "Moderate η"
    elseif st == 3
        nme = "High η"
    else
        nme = "Mixed η"
    end
    # Preallocate vector to store η values
    ηs = []
    # Preallocate vector of dissipation rates
    dsp = zeros(rps)
    # Now loop over repeats
    for i = 1:rps
        # First check that files exists
        pfile = "Data/Type$(st)/ParasType$(st)Run$(i).jld"
        if ~isfile(pfile)
            error("Run $(i) is missing a parameter file")
        end
        ofile =  "Data/Type$(st)/OutputType$(st)Run$(i).jld"
        if ~isfile(ofile)
            error("Run $(i) is missing an output file")
        end
        # Then load in data
        ps = load(pfile,"ps")
        out = load(ofile,"out")
        # Save number of surviving species
        srv[i] = ps.N
        # Store all η values
        for j = 1:ps.N
            ηs = cat(ηs,ps.mics[j].η,dims=1)
        end
        dsp[i] = dissipation(ps,out)
    end
    # Now move onto plotting
    println("Data read in")
    pyplot(dpi=200)
    # Save number of surviving species as a histogram
    histogram(srv,label="")
    plot!(title=nme,xlabel="Number of surving strains")
    savefig("Output/Type$(st)SvHist.png")
    histogram(ηs,label="")
    plot!(title=nme,xlabel=L"\eta")
    savefig("Output/Type$(st)ηHist.png")
    histogram(dsp,label="")
    plot!(title=nme,xlabel="Dissipation rate")
    savefig("Output/Type$(st)DispHist.png")
    return(nothing)
end

@time interpret()
