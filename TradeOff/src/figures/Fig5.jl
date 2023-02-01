# Script to plot figure 5
using TradeOff
using Plots
using JLD
using ColorSchemes
using Plots.PlotMeasures
using LaTeXStrings

# Function to auto generate labels
function find_label(sim_type::Int64, class::Int64)
    # Assign label based on simulation type
    if sim_type == 1
        # Check whether we are comparing free energy or maintenance
        if class == 1
            lb = "high free-energy"
        else
            lb = "low maintenance"
        end
    elseif sim_type == 2
        if class == 1
            lb = "low free-energy"
        else
            lb = "low maintenance"
        end
    elseif sim_type == 3
        if class == 1
            lb = "high free-energy"
        else
            lb = "high maintenance"
        end
    else
        if class == 1
            lb = "low free-energy"
        else
            lb = "high maintenance"
        end
    end
    return (lb)
end

function figure5(ims::Int64)
    println("Compiled")
    # Set default plotting options
    default(dpi = 200)
    # Load in colour scheme
    a = ColorSchemes.Dark2_4.colors
    # Make box for the subplots
    box = (1, bbox(0.5, 0.65, 0.4, 0.3, :bottom, :left))
    # Define reused latex strings
    e6 = L"10^6"
    # Make plot objects
    p1 = plot(xlabel = "Time (s)",
              xlim = (-Inf, 5e7),
              ylim = (0.35, 0.55),
              legend = :topleft,
              title = "Variation with free energy",
              ylabel = "Maximum ribosome fraction factor ($(L"\omega"))")
    # Add inset box to plot other trade-off into
    plot!(p1,
          xlabel = "Time ($(e6) s)",
          ylabel = L"\phi_R",
          inset_subplots = box,
          subplot = 2,
          grid = false,
          legend = false)
    plot!(p1, subplot = 2, xlim = (-Inf, 10.0), ylim = (0.05, 0.3), legend = false)
    # Now make second plot
    p2 = plot(xlabel = "Time (s)",
              xlim = (-Inf, 5e7),
              ylim = (0.4, 0.55),
              legend = :topleft,
              title = "Variation with maintenance cost",
              ylabel = "Maximum ribosome fraction factor ($(L"\omega"))")
    # Add the same inset box, to plot the other trade-off into
    plot!(p2,
          xlabel = "Time ($(e6) s)",
          ylabel = L"\phi_R",
          inset_subplots = box,
          subplot = 2,
          grid = false)
    plot!(p2, subplot = 2, xlim = (-Inf, 10.0), ylim = (0.05, 0.3), legend = false)
    p3 = plot(xlabel = "Time (s)", xlim = (-Inf, 5e7), legend = :topleft)
    p4 = plot(xlabel = "Time (s)", xlim = (-Inf, 5e7), legend = :topleft)
    # Loop over the 3 conditions
    for i in 1:3
        # Extract other simulation parameters from the function
        Np, Nt, M, d, μrange = sim_paras(i)
        # Find file name to load in
        tfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/RunStats$(ims)Ims.jld"
        # Check it actually exists
        if ~isfile(tfile)
            error("missing stats file for $(ims) immigrations simulations")
        end
        # Now load out the times, and number of trajectories
        times = load(tfile, "times")
        no_via = load(tfile, "no_via")
        # Load in averages
        mn_via_ω_bw = load(tfile, "mn_via_ω_bw")
        mn_via_ϕR = load(tfile, "mn_via_ϕR")
        mn_via_a = load(tfile, "mn_via_a")
        # Load in standard deviations
        sd_via_ω_bw = load(tfile, "sd_via_ω_bw")
        sd_via_ϕR = load(tfile, "sd_via_ϕR")
        sd_via_a = load(tfile, "sd_via_a")
        # Calculate relevant standard errors
        se_via_ω_bw = sd_via_ω_bw ./ sqrt.(no_via)
        se_via_ϕR = sd_via_ϕR ./ sqrt.(no_via)
        se_via_a = sd_via_a ./ sqrt.(no_via)
        # Plot the data to the relevant plot objects
        if i == 1 || i == 2
            # Find appropriate label
            lb = find_label(i, 1)
            plot!(p1, times, mn_via_ω_bw, ribbon = se_via_ω_bw, label = lb, color = a[i])
            plot!(p1, times / 1e6, mn_via_ϕR, ribbon = se_via_ϕR, color = a[i], subplot = 2)
            plot!(p3, times, mn_via_a, ribbon = se_via_a, label = lb, color = a[i])
        end
        if i == 1 || i == 3
            # Find appropriate label
            lb = find_label(i, 2)
            plot!(p2, times, mn_via_ω_bw, ribbon = se_via_ω_bw, label = lb, color = a[i])
            plot!(p2, times / 1e6, mn_via_ϕR, ribbon = se_via_ϕR, color = a[i], subplot = 2)
            plot!(p4, times, mn_via_a, ribbon = se_via_a, label = lb, color = a[i])
        end
    end
    # Check if directory exists and if not make it
    if ~isdir("Output/Fig5")
        mkdir("Output/Fig5")
    end
    # Save figures to this directory
    savefig(p1, "Output/Fig5/omega_with_DG.png")
    savefig(p2, "Output/Fig5/omega_with_d.png")
    savefig(p3, "Output/Fig5/a_with_DG.png")
    savefig(p4, "Output/Fig5/a_with_d.png")
    # Add annotations
    px, py = annpos([0.0; 5e7], [0.35; 0.55], 0.075, 0.0)
    annotate!(p1, px, py, text("A", 17, :black))
    px, py = annpos([0.0; 5e7], [0.4; 0.55], 0.075, 0.0)
    annotate!(p2, px, py, text("B", 17, :black))
    # Plot all graphs as a single figure
    pt = plot(p1, p2, layout = (2, 1), size = (600, 800), margin = 5.0mm)
    savefig(pt, "Output/Fig5/figure5.png")
    return (nothing)
end

@time figure5(500)
