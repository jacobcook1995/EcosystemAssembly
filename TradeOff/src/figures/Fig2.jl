# Script to plot elements needed for figure 1
using TradeOff
using Plots
using JLD2

# Function to plot the efficiency with changing ribosome fraction
function early_immigration_dyns(rN::Int64, ims::Int64, sim_type::Int64)
    println("Compiled")
    # Extract other simulation parameters from the function
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    data_dir = joinpath(
        pwd(), "Output", "$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)")
    # Read in appropriate files
    parameter_file = joinpath(data_dir, "Paras$(ims)Ims.jld")
    if ~isfile(parameter_file)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    output_file = joinpath(data_dir, "Run$(rN)Data$(ims)Ims.jld")
    if ~isfile(output_file)
        error("$(ims) immigrations run $(rN) is missing an output file")
    end
    # Read in relevant data
    ps = load(parameter_file, "ps")
    traj = load(output_file, "traj")
    T = load(output_file, "T")
    micd = load(output_file, "micd")
    its = load(output_file, "its")
    println("Data read in")
    # Find C from a function
    C = merge_data(ps, traj, T, micd, its)
    println("Data merged")
    # Define output directory and if necessary make it
    outdir = joinpath(pwd(), "Output", "Fig2")
    mkpath(outdir)
    # Find total number of strains
    totN = length(micd)
    # Set maximum time to plot to
    Tmax = 5e5
    # Set default plotting options
    default(dpi = 200)
    # Plot all the populations
    p1 = plot(yaxis = :log10,
        ylabel = "Population (# cells)",
        ylims = (1e-5, Inf),
        xlabel = "Time (s)")
    for i in 1:totN
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:, i] .> 0) .& (T .<= Tmax)
        plot!(p1, T[inds], C[inds, i], label = "")
    end
    savefig(p1, joinpath(outdir, "all_pops.png"))
    return (nothing)
end

@time early_immigration_dyns(1, 500, 1)
