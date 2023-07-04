# Script to plot figure 7 (which shows situation without immigration)
using TradeOff
using Plots
using JLD
using ColorSchemes
using LaTeXStrings
using Plots.PlotMeasures

function sub_prob(M::Int64, R::Int64, subs::Float64)
    # Check that final waste product hasn't been generated
    if subs >= M - 1
        no_sub = M - 1
        # Also check that it hasn't somehow gone negative
    elseif subs <= 0.0
        no_sub = 0.0
    else
        no_sub = subs
    end
    # Calculate probability
    P = 1 - (1 - (no_sub) / (M - 1))^R
    return (P)
end

# Function to plot the 7th figure (e.g. development without immigration)
function figure7()
    println("Compiled")
    # No immigration situation (aka sim type 5)
    ims = 0
    sim_type = 5
    tk = "NoImm"
    # Extract other simulation parameters from the function
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Find file name to load in
    sfile = "Output/$(tk)$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/RunStats$(ims)Ims.jld"
    # Check it actually exists
    if ~isfile(sfile)
        error("missing stats file for $(ims) immigrations simulations")
    end
    # Now load out the times, and number of trajectories
    times = load(sfile, "times")
    no_sims = load(sfile, "no_simulations")
    no_via = load(sfile, "no_viable_simulations")
    # Load in averages
    mn_sbs = load(sfile, "mean_no_substrates")
    mn_via_R = load(sfile, "mean_viable_species_per_reac_class")
    mn_KS_R = load(sfile, "mean_average_KS_per_reac_class")
    # Load in standard deviations
    sd_sbs = load(sfile, "sd_no_substrates")
    sd_via_R = load(sfile, "sd_viable_species_per_reac_class")
    sd_KS_R = load(sfile, "sd_average_KS_per_reac_class")
    # Preallocate standard errors
    se_via_R = zeros(size(sd_via_R))
    se_KS_R = zeros(size(sd_KS_R))
    # Calculate standard errors from this
    se_sbs = sd_sbs ./ sqrt.(no_sims)
    # Calculation (slightly) different in the viable case
    for i in axes(sd_via_R, 1)
        se_via_R[i, :] = sd_via_R[i, :] ./ sqrt.(no_via)
        se_KS_R[i, :] = sd_KS_R[i, :] ./ sqrt.(no_via)
    end
    # Preallocate probabilities
    mn_Ps = zeros(size(mn_via_R))
    up_Ps = zeros(size(mn_via_R))
    dw_Ps = zeros(size(mn_via_R))
    # Loop over, calculating the probability at each step
    for i in axes(mn_Ps, 1)
        for j in axes(mn_Ps, 2)
            # Calculate mean probability
            mn_Ps[i, j] = sub_prob(M, i, mn_sbs[j])
            # And also upper and lower bound
            up_Ps[i, j] = sub_prob(M, i, mn_sbs[j] + se_sbs[j]) .- mn_Ps[i, j]
            dw_Ps[i, j] = mn_Ps[i, j] .- sub_prob(M, i, mn_sbs[j] - se_sbs[j])
        end
    end
    println("Data read in")
    # Check if directory exists and if not make it
    if ~isdir("Output/Fig7")
        mkdir("Output/Fig7")
    end
    # Set default plotting options
    default(dpi = 200)
    # Load in colour scheme
    a = ColorSchemes.sunset.colors
    # Plot basic trade-off first
    p1 = plot(xlabel = "Time (s)", ylabel = "Number of species", xlim = (-Inf, 2.5e6))
    plot!(p1,
          title = "Number of reactions with time",
          legend = :bottomright,
          ylim = (0.0, 5.0))
    plot!(p1, times, mn_via_R[1, :], ribbon = se_via_R[1, :], label = "R=1", color = a[1])
    plot!(p1, times, mn_via_R[3, :], ribbon = se_via_R[3, :], label = "R=3", color = a[2])
    plot!(p1, times, mn_via_R[5, :], ribbon = se_via_R[5, :], label = "R=5", color = a[3])
    plot!(p1, times, mn_via_R[7, :], ribbon = se_via_R[7, :], label = "R=7", color = a[4])
    # Add annotation
    px, py = annpos([0.0; 2.5e6], [0.0; 5.0], 0.075, 0.05)
    annotate!(p1, px, py, text("A", 17, :black))
    savefig(p1, "Output/Fig7/AvViaReacsTime.png")
    # Now do probability plot
    p2 = plot(xlabel = "Time (s)",
              ylabel = "Probability of no usable substrate",
              xlim = (-Inf, 2.5e6),
              title = "Chance of species finding no usable substrates",
              legend = false,
              times,
              1 .- mn_Ps[1, :],
              ribbon = (up_Ps[1, :], dw_Ps[1, :]),
              label = "R=1",
              color = a[1])
    plot!(p2,
          times,
          1 .- mn_Ps[3, :],
          ribbon = (up_Ps[3, :], dw_Ps[3, :]),
          label = "R=3",
          color = a[2])
    plot!(p2,
          times,
          1 .- mn_Ps[5, :],
          ribbon = (up_Ps[5, :], dw_Ps[5, :]),
          label = "R=5",
          color = a[3])
    plot!(p2,
          times,
          1 .- mn_Ps[7, :],
          ribbon = (up_Ps[7, :], dw_Ps[7, :]),
          label = "R=7",
          color = a[4])
    # Define box for inset here
    box = (1, bbox(0.5, 0.5, 0.45, 0.3375, :bottom, :left))
    Ks = L"K_S"
    e7 = L"10^7"
    em3 = L"10^{-3}"
    plot!(p2,
          times / 1e7,
          mn_KS_R[1, :] * 1000.0,
          ribbon = se_KS_R[1, :] * 1000.0,
          label = "R=1",
          color = a[1],
          inset_subplots = box,
          subplot = 2)
    plot!(p2,
          times / 1e7,
          mn_KS_R[3, :] * 1000.0,
          ribbon = se_KS_R[3, :] * 1000.0,
          label = "R=3",
          color = a[2],
          subplot = 2)
    plot!(p2,
          times / 1e7,
          mn_KS_R[5, :] * 1000.0,
          ribbon = se_KS_R[5, :] * 1000.0,
          label = "R=5",
          color = a[3],
          subplot = 2)
    plot!(p2,
          times / 1e7,
          mn_KS_R[7, :] * 1000.0,
          ribbon = se_KS_R[7, :] * 1000.0,
          label = "R=7",
          color = a[4],
          subplot = 2)
    plot!(p2,
          xlim = (-Inf, 0.25),
          guidefontsize = 9,
          grid = false,
          legend = false,
          subplot = 2)
    plot!(p2, xlabel = "Time ($(e7) s)", ylabel = "$(Ks) ($em3)", subplot = 2)
    # Add annotation
    px, py = annpos([0.0; 2.5e6], [0.0; 1.1], 0.075, 0.05)
    annotate!(p2, px, py, text("B", 17, :black))
    savefig(p2, "Output/Fig7/ProbSubTime.png")
    # Plot all graphs as a single figure
    pt = plot(p1, p2, layout = (2, 1), size = (600, 800), margin = 5.0mm)
    savefig(pt, "Output/Fig7/figure7.png")
    return (nothing)
end

@time figure7()
