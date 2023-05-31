# Script to plot figure 3
using TradeOff
using Plots
using JLD
using ColorSchemes
using LaTeXStrings
using Plots.PlotMeasures

function figure3(ims::Int64, sim_type::Int64)
    println("Compiled")
    # Extract other simulation parameters from the function
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Read in appropriate files
    lfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/SurvTimes$(ims)Ims.jld"
    if ~isfile(lfile)
        error("$(ims) immigrations simulation is missing a long term survivors file")
    end
    # Read in relevant data
    sTs = load(lfile, "sTs")
    # Now look at snapshot data
    sfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/SnapDataStats$(ims)Ims.jld"
    if ~isfile(sfile)
        error("$(ims) immigrations simulation is missing an snapshot stats file")
    end
    # Save growth probabilities
    gp = load(sfile, "gp")
    s_times = load(sfile, "times")
    # Find file name to load in
    tfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/RunStats$(ims)Ims.jld"
    # Check it actually exists
    if ~isfile(tfile)
        error("missing stats file for $(ims) immigrations simulations")
    end
    # Now load out the times, and number of trajectories
    t_times = load(tfile, "times")
    no_simulations = load(tfile, "no_simulations")
    # Load in final ϕR values
    all_final_ϕRs = load(tfile, "all_final_ϕRs")
    # Load in averages
    mean_surviving_species = load(tfile, "mean_surviving_species")
    mean_total_population = load(tfile, "mean_total_population")
    mean_shannon_diversity = load(tfile, "mean_shannon_diversity")
    # Load in standard deviations
    sd_surviving_species = load(tfile, "sd_surviving_species")
    sd_total_population = load(tfile, "sd_total_population")
    sd_shannon_diversity = load(tfile, "sd_shannon_diversity")
    # Calculate standard errors
    se_surviving_species = sd_surviving_species ./ sqrt.(no_simulations)
    se_total_population = sd_total_population ./ sqrt.(no_simulations)
    se_shannon_diversity = sd_shannon_diversity ./ sqrt.(no_simulations)
    println("Data read in")
    # Check if directory exists and if not make it
    if ~isdir("Output/Fig3")
        mkdir("Output/Fig3")
    end
    # Set default plotting options
    default(dpi = 200)
    # Define latex commands
    e7 = L"10^7"
    e13 = L"10^{13}"
    # Load in colour scheme
    a = ColorSchemes.tab10.colors
    # Make first plot
    p1 = plot(ylabel = "Total population ($(e13) cells)",
              xlim = (-Inf, 5e7),
              legend = :topleft,
              t_times,
              mean_total_population / 1e13,
              ribbon = se_total_population / 1e13,
              label = "Population",
              color = a[1],
              ylim = (-Inf, 4.0))
    # Define box for inset here
    box = (1, bbox(0.4, 0.2, 0.4, 0.3, :bottom, :left))
    # Add histogram in as an insert
    histogram!(p1,
               all_final_ϕRs,
               nbins = 100,
               color = a[2],
               label = "",
               inset_subplots = box,
               xlabel = "Final ribosome fraction ($(L"\phi_R"))",
               ylabel = "Number of survivors",
               subplot = 2)
    plot!(p1, guidefontsize = 9, grid = false, subplot = 2)
    # Twin the x-axis
    pt = twinx(p1)
    # Ensure same limits are used
    plot!(pt, xlim = (-Inf, 5e7), ylim = (-Inf, 2.9), ylabel = "Shannon diversity")
    # Then plot the Shannon diversity
    p1 = plot!(pt,
               t_times,
               mean_shannon_diversity,
               ribbon = se_shannon_diversity,
               label = "Diversity",
               color = a[3],
               legend = :topright)
    # Add annotation
    px, py = annpos([0.0; 5e7], [0.0; 3.9], 0.05, 0.0)
    annotate!(p1, px, py, text("A", 17, :black))
    savefig(p1, "Output/Fig3/SumStats.png")
    # Now make the second plot
    p2 = plot(xlabel = "Time (s)",
              ylabel = "Number of species",
              xlim = (-Inf, 5e7),
              ylim = (-Inf, 30.0),
              t_times,
              mean_surviving_species,
              ribbon = se_surviving_species,
              label = "Species",
              color = a[4],
              legend = :topleft)
    # Add annotation
    px, py = annpos([0.0; 5e7], [0.0; 30.0], 0.05, 0.0)
    annotate!(p2, px, py, text("B", 17, :black))
    # Define box for inset here
    box = (1, bbox(0.4, 0.4, 0.4, 0.3, :bottom, :left))
    # Add histogram in as an insert
    histogram!(p2,
               sTs / 1e7,
               nbins = 100,
               color = a[5],
               label = "",
               inset_subplots = box,
               xlabel = "Time of immigration ($(e7) s)",
               ylabel = "Number of survivors",
               subplot = 2)
    plot!(p2, guidefontsize = 9, grid = false, subplot = 2)
    # Twin the x-axis
    pt = twinx(p2)
    # Ensure same limits are used
    plot!(pt, xlim = (-Inf, 5e7), ylabel = "Proportion of immigrants that can grow")
    # Then plot the proportion feasible
    p2 = plot!(pt, s_times, gp, label = "Invasibility", color = a[6], legend = :topright)
    savefig(p2, "Output/Fig3/Invasibility.png")
    # Now want to make a plot incorporating both previous plots
    pt = plot(p1, p2, layout = (2, 1), size = (750, 800), margin = 5.0mm,
              rightmargin = 20mm)
    savefig(pt, "Output/Fig3/figure3.png")
    return (nothing)
end

@time figure3(500, 1)
