# A script to read in parameter files and test why
using Assembly
using Plots
using JLD
using LaTeXStrings
import PyPlot

function test()
    # Check that sufficient arguments have been provided
    if length(ARGS) < 2
        error("need to specify number of reactions and run number for file name")
    end
    # Preallocate the variables I want to extract from the input
    R = 0
    rN = 0
    # Check that all arguments can be converted to integers
    try
        R = parse(Int64, ARGS[1])
        rN = parse(Int64, ARGS[2])
    catch e
        error("all three inputs must be integer")
    end
    # Check that simulation type is valid
    if R < 1
        error("invalid number of reactions each strain must have more than 1 reaction")
    end
    # Check that number of strains is greater than 0
    if rN < 1
        error("need to do at least 1 simulation")
    end
    # Now check if file exists
    pfile = "Data/Temp/ParasType$(R)Run$(rN).jld"
    if ~isfile(pfile)
        error("the file you are looking for doesn't exist")
    end
    # Then load in parameters
    ps = load(pfile, "ps")
    # Now start actual script
    println("Compiled and input read in!")
    # Fairly arbitrary initial conditions
    pop = ones(ps.N)
    conc = zeros(ps.M)
    as = 1e5 * ones(ps.N)
    ϕs = 0.1 * ones(ps.N)
    # Choose a really short time span
    Tmax = 5e5
    # Then run the simulation
    C, T = test_full_simulate(ps, Tmax, pop, conc, as, ϕs)
    # Give maximum and minimum for each variable type
    println("Max population: $(maximum(C[end,1:ps.N]))")
    println("Min population: $(minimum(C[end,1:ps.N]))")
    println("Max concentration: $(maximum(C[end,(ps.N+1):(ps.N+ps.M)]))")
    println("Min concentration: $(minimum(C[end,(ps.N+1):(ps.N+ps.M)]))")
    println("Max energy: $(maximum(C[end,(ps.N+ps.M+1):(2*ps.N+ps.M)]))")
    println("Min energy: $(minimum(C[end,(ps.N+ps.M+1):(2*ps.N+ps.M)]))")
    println("Max fraction: $(maximum(C[end,(2*ps.N+ps.M+1):end]))")
    println("Min fraction: $(minimum(C[end,(2*ps.N+ps.M+1):end]))")
    # Doing pretty naive plotting to begin with
    # Setup population plot
    p1 = plot(xlabel = "Time", ylabel = "Population", yaxis = :log10)
    for i in 1:(ps.N)
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:, i] .> 0)
        plot!(p1, T[inds], C[inds, i], label = "")
    end
    savefig(p1, "Output/PopvsTime.png")
    plot(T, C[:, (ps.N + 1):(ps.N + ps.M)], xlabel = "Time", label = "",
         ylabel = "Concentration")
    savefig("Output/MetabolitevsTime.png")
    plot(T, C[:, (ps.N + ps.M + 1):(2 * ps.N + ps.M)], xlabel = "Time", label = "",
         ylabel = "Cell energy conc")
    savefig("Output/EnergyvsTime.png")
    plot(T, C[:, (2 * ps.N + ps.M + 1):end], xlabel = "Time", label = "",
         ylabel = L"\phi_R")
    savefig("Output/FractionvsTime.png")
    # Sensible plotting will require removal of dead species
    # Metabolites comparatively easy to plot though
    # Questionable if I want to show all the fractions on the same plot, they seem "internal"
    return (nothing)
end

@time test()
