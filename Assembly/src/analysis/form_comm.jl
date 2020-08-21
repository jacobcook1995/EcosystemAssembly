# A script to form stable communities using the full model that then saves them for later use
using Assembly
using Plots
using LaTeXStrings
using JLD
using SymPy
import PyPlot

# # function to test that the new stuff I'm writing actually works
# # Keeping this for the moment as the plotting stuff might be useful
# function test()
#     println("Compiled!")
#     # Assume that half saturation occurs at a quarter κ/δ
#     KS = (1/4)*5.5e-3
#     # From wikipedia an average enzyme has a k to KS ratio of 10^5 M^-1 s^-1
#     # This would give us a k of 137.5, sensible to assume an above average rate
#     # Though should be reduced by the fact we include uptake as well as metabolism
#     # Choosing k = 500 means we match the maximum glucose uptake rate seen in Natarajan et al (2000)
#     # of 3*10^7 molecules per second.
#     # The above is a sensible argument but 1.0 gives a more reasonable ATP concentration.
#     kc = 10.0
#     # The reversibility factor remains the same as previously
#     kr = 10.0
#     # Assume microbes have 2 reactions each
#     mR = 2.0
#     sdR = 0.0
#     # Case of 8 metabolites and 1 strain
#     N = 3
#     M = 8
#     # Use formula to find how many reactions this implies
#     O = 2*M - 3
#     good = false
#     # Make parameter set
#     ps = 0
#     while good == false
#         ps = initialise(N,M,O,mR,sdR,kc,KS,kr)
#         # This step is to ensure that the first metabolite can be broken down
#         for i = 1:N
#             if any((ps.reacs[ps.mics[i].Reacs].↦:Rct) .== 1)
#                 good = true
#             end
#         end
#     end
#     # Set time long enough to see dynamics
#     Tmax = 10000000.0
#     # Fairly arbitary inital conditions
#     pop = ones(N)
#     conc = zeros(M)
#     as = 1e5*ones(N)
#     ϕs = 0.1*ones(N)
#     C, T = full_simulate(ps,Tmax,pop,conc,as,ϕs)
#     pyplot(dpi=200)
#     # Calculate growth rates
#     λ = zeros(length(C[:,1]),N)
#     for i = 1:length(C[:,1])
#         for j = 1:N
#             λ[i,j] = λs(C[i,N+M+j],C[i,2*N+M+j],ps.mics[j])
#         end
#     end
#     # Then plot growth rates
#     plot(T,λ,xlabel="Time",label="",ylabel=L"\lambda")
#     savefig("Output/GrowthvsTime.png")
#     # Calculate reaction inhibition
#     # Find surviving microbes
#     ms = []
#     for i = 1:N
#         if C[end,i] > 1e-8
#             ms = cat(ms,i,dims=1)
#         end
#     end
#     # Make array of bools to store results
#     # To generalise this I would have to figure out how to make
#     # an array with rows of differing lengths
#     rs = fill(false,(length(ms),ps.mics[1].R))
#     # Then for these microbes find viable reactions
#     for i = 1:length(ms)
#         # Find corresponding microbes
#         mic = ps.mics[ms[i]]
#         for j = 1:mic.R
#             # find reaction
#             r = ps.reacs[mic.Reacs[j]]
#             if C[end,N+r.Rct] > 1e-8 && C[end,N+r.Prd] > 1e-8
#                 rs[i,j] = true
#             end
#         end
#     end
#     # Preallocate array of θs, this also will take some work to generalise
#     θs = zeros(length(C[:,1]),size(rs,1)*size(rs,2))
#     for i = 1:length(ms)
#         for j = 1:ps.mics[1].R
#             # Find appropriate reaction
#             r = ps.reacs[ps.mics[ms[i]].Reacs[j]]
#             for k = 1:length(C[:,1])
#                 θs[k,(i-1)*ps.mics[1].R+j] = θ_smooth(C[k,N+r.Rct],C[k,N+r.Prd],ps.T,ps.mics[ms[i]].η[j],r.ΔG0)
#             end
#         end
#     end
#     # Then plot inhibitions
#     plot(T,θs,xlabel="Time",label="",ylabel=L"\theta")
#     savefig("Output/InhibvsTime.png")
#     return(nothing)
# end

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
    Tmax = 1e7
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

# function to calculate the dissipation for an assembled ecosystem
function dissipation(ps::FullParameters,out::Array{Float64,1})
    # check that parameter set is sensible given the output
    if length(out) != ps.M + 3*ps.N
        error("parameter set doesn't match output")
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
            # Check if reaction actually occurs
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
    println("Compiled!")
    # Preallocate memory to store the number of survivors in
    srv = zeros(rps)
    # Preallocate vector to store percentage of free energy dissipated (under standard conditions)
    pd = []
    # Preallocate vector of dissipation rates
    dsp = zeros(rps)
    # Now loop over repeats
    for i = 1:rps
        # First check that files exists
        pfile = "Data/Type$(R)/ParasType$(R)Run$(i).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        ofile = "Data/Type$(R)/OutputType$(R)Run$(i).jld"
        if ~isfile(ofile)
            error("run $(i) is missing an output file")
        end
        efile = "Data/Type$(R)/ExtinctType$(R)Run$(i).jld"
        if ~isfile(ofile)
            error("run $(i) is missing an extinct file")
        end
        # Then load in data
        ps = load(pfile,"ps")
        out = load(ofile,"out")
        # Save number of surviving species
        srv[i] = ps.N
        # Store all percentage dissipations
        for j = 1:ps.N
            # Find vector of η values
            ηs = ps.mics[j].η
            # make temporary vector for percentages
            pdt = zeros(length(ηs))
            # loop over η values
            for k = 1:length(ηs)
                # Find standard Gibbs free energy
                dG = ps.reacs[ps.mics[j].Reacs[k]].ΔG0
                # Calculate percentage dissipated (under standard conditions)
                pdt[k] = -dG/(ηs[k]*ΔGATP)
            end
            pd = cat(pd,pdt,dims=1)
        end
        # Remove any negative metabolite concentrations
        for i = (ps.N+1):(ps.N+ps.M)
            if out[i] < 0.0
                out[i] = 0.0
            end
        end
        # Now its safe to calculate the dissipation
        dsp[i] = dissipation(ps,out)
        # Going to do this only for the final step for now
        if i == rps
            # Preallocate vector of forces
            F = Array{Sym,1}(undef,3*ps.N+ps.M)
            # Find forces using function
            F = Force(ps,F)
            # Use final simulation results to find local forces
            f = nForce(F,out,ps)
            # Preallocate necessary data structures
            dx = zeros(length(f))
            rate = zeros(ps.N,ps.O)
            t = 0.0 # No explict time dependance so this doesn't matter
            dx = test_dynamics!(dx,out,ps,rate,t)
            println(f)
            println(dx)
            println(f./dx)
            # # Check maximum force
            # stabN = (maximum(abs.(f[1:ps.N])) <= 2.0)
            # println("Populations stable = $(stabN)")
            # println(abs.(f[1:ps.N]))
            # # No seems to work now
            # stabM = (maximum(abs.(f[ps.N+1:ps.N+ps.M])) <= 1.0e-5)
            # println("Concentrations stable = $(stabM)")
            # println(abs.(f[ps.N+1:ps.N+ps.M]))
            # staba = (maximum(abs.(f[(ps.N+ps.M+1):(2*ps.N+ps.M)])) <= 1.0e-5)
            # println("Energies stable = $(staba)")
            # println(abs.(f[(ps.N+ps.M+1):(2*ps.N+ps.M)]))
            # stabϕ = (maximum(abs.(f[(2*ps.N+ps.M+1):end])) <= 1.0e-5)
            # println("Fractions stable = $(stabϕ)")
            # println(abs.(f[(2*ps.N+ps.M+1):end]))
            # # Find corresponding data
            # T = load(ofile,"T")
            # C = load(ofile,"C")
            # # Useful to also have the extinction data available
            # ded = load(efile,"ded")
            # # Setup population plot
            # pyplot(dpi=200)
            # p1 = plot(xlabel="Time",ylabel="Population",yaxis=:log10)
            # # Find number of ignored microbes
            # N = ps.N + length(ded)
            # for i = 1:N
            #     # Find and eliminate zeros so that they can be plotted on a log plot
            #     inds = (C[:,i] .> 0)
            #     plot!(p1,T[inds],C[inds,i],label="")
            # end
            # savefig(p1,"Output/PopvsTime.png")
            # # Now plot the metabolites
            # plot(T,C[:,(N+1):(N+ps.M)],xlabel="Time",label="",ylabel="Concentration")
            # savefig("Output/MetabolitevsTime.png")
            # # Now plot the energy concentrations
            # plot(T,C[:,(N+ps.M+1):(2*N+ps.M)],xlabel="Time",label="",ylabel="Cell energy conc")
            # savefig("Output/EnergyvsTime.png")
            # plot(T,C[:,2*N+ps.M+1:end],xlabel="Time",label="",ylabel=L"\phi_R")
            # savefig("Output/FractionvsTime.png")
        end
    end
    # Now move onto plotting
    println("Data read in")
    pyplot(dpi=200)
    # Save number of surviving species as a histogram
    histogram(srv,label="")
    plot!(title="$(R) reactions",xlabel="Number of surving strains")
    savefig("Output/Type$(R)SvHist.png")
    histogram(pd,label="")
    plot!(title="$(R) reactions",xlabel="Percentage of free energy dissipated")
    savefig("Output/Type$(R)ηHist.png")
    histogram(dsp,label="")
    plot!(title="$(R) reactions",xlabel="Dissipation rate")
    savefig("Output/Type$(R)DispHist.png")
    return(nothing)
end

@time interpret()
