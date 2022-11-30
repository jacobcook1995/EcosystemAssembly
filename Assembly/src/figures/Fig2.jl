# Script to construct figure 2
using Assembly
using Plots
using JLD
using StatsBase
using KernelDensity
using StatsPlots
using Plots.PlotMeasures
using LsqFit
import PyPlot

# function to calculate the dissipation for an assembled ecosystem
function dissipation(ps::FullParameters, ms::Array{MicrobeP, 1}, out::Array{Float64, 1})
    # Set all elements of out less than zero to zero
    out[out .< 0.0] .= 0.0
    # Define number of strains
    N = length(ms)
    # check that parameter set is sensible given the output
    if length(out) != ps.M + 3 * N
        error("parameter set doesn't match output")
    end
    # Set dissipation to zero
    dsp = 0
    # Loop over number of strains
    for i in 1:N
        # Isolate this strain
        mic = ms[i]
        # Loop over reactions of this strain
        for j in 1:(mic.R)
            # Find appropriate reaction
            r = ps.reacs[mic.Reacs[j]]
            # If there's no product Gibbs free energy becomes infinite. Justified to ignore
            # this as if product hasn't built up reaction can't be happening to a significant degree
            if out[N + r.Prd] != 0.0
                # Find amount of energy that this reaction dissipates
                Fd = -(r.ΔG0 + Rgas * ps.T * log(out[N + r.Prd] / out[N + r.Rct]) +
                       mic.η[j] * ΔGATP)
                # Find amount of enzyme E
                E = Eα(out[2 * N + ps.M + i], mic, j)
                # Then find the rate that this reaction proceeds at
                q = qs(out[N + r.Rct], out[N + r.Prd], E, j, mic, ps.T, r)
                # Check if reaction actually occurs
                if q != 0.0
                    dsp += q * Fd * out[i]
                end
            end
        end
    end
    # Convert from molecule units to moles
    dsp /= NA
    return (dsp)
end

# Function to count the number of spikes in the entropy production trace
function ent_sp_cnt(ps::FullParameters, C::Array{Float64, 2}, T::Array{Float64, 1},
                    Ns::Int64,
                    out::Array{Float64, 1}, inf_out::Array{Float64, 1},
                    ded::Array{MicrobeP, 1})
    # Make new vector of microbes
    ms = Array{MicrobeP, 1}(undef, Ns)
    # Setup counter
    cnt = 0
    # Loop over all microbes
    for i in 1:Ns
        # Check if it is a survivor
        if C[end, i] != 0.0 && C[end, i] ∈ out
            # If it is find and save it
            ind = findfirst(x -> x == C[end, i], out)
            ms[i] = ps.mics[ind]
        else
            # Update counter
            cnt += 1
            # Use next element from ded vector
            ms[i] = ded[cnt]
        end
    end
    # Find final time
    Tmax = T[end]
    # container to store entropy production
    ep = zeros(length(T))
    # Calculate entropy production at each step
    for i in 1:length(T)
        # Calculate entropy production at each step
        ep[i] = dissipation(ps, ms, C[i, :])
    end
    # container to store changes in entropy productions
    dep = zeros(length(T))
    # Find change between one point and the previous one
    for i in 1:length(T)
        if i > 1
            dep[i] = ep[i] - ep[i - 1]
        else
            dep[i] = ep[i]
        end
    end
    # Find first peak, i.e. first time entropy production decays
    pind = findfirst(x -> x < 0.0, dep) - 1
    # Store in a vector
    pinds = [pind]
    # Counter for peaks
    pc = 1
    # Has the end of the trajectory been reached
    nd = false
    # Set initial offset
    off = pind
    # Now identify every peak
    while nd == false
        # Find first point where entropy production increases again
        try
            tind = findfirst(x -> x > 0.0, dep[(off + 1):end]) - 1
            # Add trough to offset
            off += tind
            # Error arises because we're at the end of the vector
        catch e
            # Set offset to maximum
            off = length(T)
            # and finish loop
            nd = true
        end
        # Then find next peak based on this
        try
            pind = findfirst(x -> x < 0.0, dep[(off + 1):end]) - 1
            # Error arises because we're at the end of the vector
        catch e
            # Set final peak
            pind = length(T) - off
            # and finish loop
            nd = true
        end
        # Check difference between peak and trough height
        df = (abs(ep[(pind + off)]) - abs(ep[(off)])) / (abs(ep[(pind + off)]))
        # If a greater than 1% change save
        if df > 0.01
            # Add new value to the vector
            pinds = cat(pinds, pind + off, dims = 1)
            # Increment counter
            pc += 1
        end
        # Finally add 'peak' to the offset, regardless of if it's a true peak or not
        off += pind
    end
    # Check number of peaks is reasonable (can be wrong due to noisy trajectories)
    if pc > 100
        # If it's not set as a NaN
        pc = NaN
    end
    return (pc)
end

function figure2(Rl::Int64, Ru::Int64, syn::Bool, Nr::Int64, Ns::Int64, en::String,
                 Tf::Float64, rps::Int64)
    println("Compiled!")
    # Set initial number of concentrations, and preallocate the final ones
    nmi = ones(rps)
    nmf = zeros(rps)
    # Preallocate number of entropy production spikes
    pcs = zeros(rps)
    # Loop over all repeats to find substrate diversification
    for i in 1:rps
        # Read in relevant files
        pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ns).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        ofile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/RedOutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ns).jld"
        if ~isfile(ofile)
            error("run $(Nr) is missing an output file")
        end
        # Load required data
        ps = load(pfile, "ps")
        inf_out = load(ofile, "inf_out")
        # Only want to count substrates
        nmf[i] = count(x -> x > 0.0, inf_out[(ps.N + 1):(ps.N + ps.M - 1)])
    end
    # Loop to count entropy production spikes
    for i in 1:rps
        # Read in specific files needed for the dynamics
        pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/ParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ns).jld"
        if ~isfile(pfile)
            error("run $(Nr) is missing a parameter file")
        end
        ofile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/OutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ns).jld"
        if ~isfile(ofile)
            error("run $(Nr) is missing an output file")
        end
        ofile2 = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/RedOutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ns).jld"
        if ~isfile(ofile)
            error("run $(Nr) is missing a 2nd output file")
        end
        efile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/ExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ns).jld"
        if ~isfile(efile)
            error("run $(Nr) is missing an extinct file")
        end
        # Read in relevant data
        ps = load(pfile, "ps")
        C = load(ofile, "C")
        T = load(ofile, "T")
        out = load(ofile, "out")
        inf_out = load(ofile2, "inf_out")
        ded = load(efile, "ded")
        # Then run peak counting function
        pcs[i] = ent_sp_cnt(ps, C, T, Ns, out, inf_out, ded)
    end
    # Read in specific files needed for the dynamics
    pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/ParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(Nr)Ns$(Ns).jld"
    if ~isfile(pfile)
        error("run $(Nr) is missing a parameter file")
    end
    ofile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/OutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(Nr)Ns$(Ns).jld"
    if ~isfile(ofile)
        error("run $(Nr) is missing an output file")
    end
    efile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/ExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(Nr)Ns$(Ns).jld"
    if ~isfile(efile)
        error("run $(Nr) is missing an extinct file")
    end
    ps = load(pfile, "ps")
    C = load(ofile, "C")
    T = load(ofile, "T")
    out = load(ofile, "out")
    ded = load(efile, "ded")
    # Make new vector of microbes
    ms = Array{MicrobeP, 1}(undef, Ns)
    # Setup counter
    cnt = 0
    # Loop over all microbes
    for i in 1:Ns
        # Check if it is a survivor
        if C[end, i] != 0.0 && C[end, i] ∈ out
            # If it is find and save it
            ind = findfirst(x -> x == C[end, i], out)
            ms[i] = ps.mics[ind]
        else
            # Update counter
            cnt += 1
            # Use next element from ded vector
            ms[i] = ded[cnt]
        end
    end
    # Find final time
    Tmax = T[end]
    # Find time to plot too
    Tend = Tmax * Tf
    # Strains to display
    ds = 25
    # Check if this is enough to show all survivors
    if ds < ps.N
        error("number of survivors higher than number of strains shown")
    end
    # Find indices of strains I want to show
    is = zeros(Int64, ds)
    for i in 1:(ps.N)
        is[i] = findfirst(x -> x == out[i], C[end, :])
    end
    # Fill in remaining non-survivors
    for i in (ps.N + 1):ds
        # Find strains not already included
        dff = setdiff(1:Ns, is[1:(i - 1)])
        is[i] = dff[1]
    end
    # Set suitable threshold for viable usage of metabolite 3
    tsh3 = 1e-4
    # Preallocate threshold times
    Tms = zeros(ps.M)
    # Find and save times where metabolites have first crossed threshold
    for i in 1:(ps.M)
        ind = findfirst(x -> x >= tsh3, C[:, Ns + i])
        # Check if threshold actually is ever crossed
        if ind != nothing
            Tms[i] = T[ind]
        else
            Tms[i] = NaN
        end
    end
    # Set threshold to be considered accumulated
    acct = 2e-3
    # Preallocate checks that they went over the threshold
    mtr = fill(false, ps.M)
    # And exhaustion times
    exT = zeros(ps.M)
    # Loop over all metabolites and check if they went over the threshold
    for i in 1:(ps.M)
        # Check if concentration goes over threshold
        mtr[i] = any(x -> x >= acct, C[:, Ns + i])
        # If it does find exhaustion time
        if mtr[i] == true
            # Find index of first going over the time
            ind = findfirst(x -> x >= acct, C[:, Ns + i])
            # Then check if it drops below again
            if any(x -> x < acct, C[(ind + 1):end, Ns + i])
                # If so find first point below the threshold
                ind2 = findfirst(x -> x < acct, C[(ind + 1):end, Ns + i])
                # Use to find exhaustion time
                exT[i] = T[ind2 + ind]
            else
                exT[i] = NaN
            end
        else
            exT[i] = NaN
        end
    end
    # Find indices of NaNs
    nans = isnan.(pcs)
    # Need to reduce data here to only include the relevant points
    xdataT = pcs[.!nans]
    ydataT = nmf[.!nans] .- 1
    # Now calculate Pearson correlation coefficient
    xbarT = sum(xdataT) / length(xdataT)
    ybarT = sum(ydataT) / length(ydataT)
    a = 0
    b = 0
    c = 0
    for i in eachindex(xdataT)
        a += (xdataT[i] - xbarT) * (ydataT[i] - ybarT)
        b += (xdataT[i] - xbarT)^2
        c += (ydataT[i] - ybarT)^2
    end
    r = a / sqrt(b * c)
    println("Correlation between peaks and substrates: $(r)")
    # Set model to fit the line to
    @. model(x, p) = p[1] + p[2] * x
    p0 = [0.0, 1.0] # Initial values
    # Fit model
    fitT = curve_fit(model, xdataT, ydataT, p0)
    # Extract values
    yintT = coef(fitT)[1]
    slopT = coef(fitT)[2]
    println("Intercept = $(yintT)")
    println("Slope = $(slopT)")
    # Set line width (plots defaults to 1)
    wdt = 2
    # Now move onto plotting
    pyplot()
    theme(:wong2, dpi = 300, guidefontsize = 16, tickfontsize = 14)
    wongc = wong2_palette()
    # Plot all the populations
    p1 = plot(yaxis = :log10, ylabel = "Population (# cells)")
    # Store max and min C values
    maxC = zeros(length(is))
    minC = zeros(length(is))
    c = 0
    for i in is
        c += 1
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:, i] .> 0) .& (T .<= Tend)
        plot!(p1, T[inds], C[inds, i], lw = wdt, label = "")
        # Store max and min C values in range
        maxC[c] = maximum(C[inds, i])
        minC[c] = minimum(C[inds, i])
    end
    # Add annotation
    px, py = annpos([0.0; Tend], [maxC; minC])
    # Different because log10 scale used
    annotate!(px, py, text("A", 20, :black))
    vline!(p1, [Tms[3]], color = :red, style = :dash, lw = wdt, label = "")
    savefig(p1, "Output/Fig2/pops.png")
    # Now plot concentrations
    p2 = plot(ylabel = "Metabolite concentration (moles)")
    # Store max and min C values
    maxC = zeros(length(is))
    minC = zeros(length(is))
    c = 0
    for i in (Ns + 1):(Ns + ps.M)
        c += 1
        # Find and eliminate points after end time
        inds = (T .<= Tend)
        # Can't switch theme but can switch palette to avoid repeated colours
        if (i - Ns) >= 20
            plot!(p2, T[inds], C[inds, i], lw = wdt, label = "", palette = :darktest)
        else
            plot!(p2, T[inds], C[inds, i], lw = wdt, label = "")
        end
        # Store max and min C values in range
        maxC[c] = maximum(C[inds, i])
        minC[c] = minimum(C[inds, i])
    end
    # Add annotation
    px, py = annpos([0.0; Tend], [maxC; minC], 0.10, 0.05)
    annotate!(px, py, text("C", 20, :black))
    vline!(p2, [Tms[3]], color = :red, style = :dash, lw = wdt, label = "")
    # Define box for inset here
    box = (1, bbox(0.65, 0.25, 0.325, 0.275, :bottom, :left))
    # Find kernel densities directly
    dei = kde(nmi, boundary = (0, ps.M - 1))
    def = kde(nmf, boundary = (0, ps.M - 1))
    # Increase final density
    def.density = 4 * def.density
    # Plot as density plots
    plot!(p2, dei, color = :black, label = "Initial", inset_subplots = box, subplot = 2)
    plot!(p2[2], def, color = :red, label = "Final", xlabel = "Number of substrates")
    plot!(p2[2], guidefontsize = 12, legendfontsize = 12, tickfontsize = 9, yaxis = false,
          grid = false)
    savefig(p2, "Output/Fig2/concs.png")
    # Now plot proteome fraction
    p3 = plot(ylabel = "Ribosome fraction")
    # Store max and min C values
    maxC = zeros(length(is))
    minC = zeros(length(is))
    c = 0
    for i in is .+ 2 * Ns .+ ps.M
        c += 1
        # Find and eliminate points after end time, remove points where strain is dead
        inds = (C[:, i - 2 * Ns - ps.M] .> 0) .& (T .<= Tend)
        plot!(p3, T[inds], C[inds, i], lw = wdt, label = "")
        # Store max and min C values in range
        maxC[c] = maximum(C[inds, i])
        minC[c] = minimum(C[inds, i])
    end
    # Add annotation
    px, py = annpos([0.0; Tend], [maxC; minC])
    annotate!(px, py, text("B", 20, :black))
    vline!(p3, [Tms[3]], color = :red, style = :dash, lw = wdt, label = "")
    savefig(p3, "Output/Fig2/fracs.png")
    # container to store entropy production
    ep = zeros(length(T))
    # Calculate entropy production at each step
    for i in 1:length(T)
        # Calculate entropy production at each step
        ep[i] = dissipation(ps, ms, C[i, :])
    end
    p4 = plot(xlabel = "Time (s)", ylabel = "Entropy production (J/K per s)")
    # Find and eliminate points after end time
    inds = (T .<= Tend)
    plot!(p4, T[inds], ep[inds], lw = wdt, label = "", ylim = (-0.01, Inf))
    # Add annotation
    px, py = annpos([0.0; Tend], ep[inds])
    annotate!(px, py, text("D", 20, :black))
    vline!(p4, [Tms[3]], color = :red, style = :dash, lw = wdt, label = "")
    for i in 1:(ps.M)
        if mtr[i] == true && ~isnan(exT[i])
            println(exT[i])
            plot!(p4, [exT[i]; exT[i]], [-0.01; 0.01], color = wongc[i], style = :solid,
                  lw = wdt, label = "")
        end
    end
    # Define box for inset here
    box = (1, bbox(0.65, 0.20, 0.325, 0.325, :bottom, :left))
    # Plot scatter points as an inset
    scatter!(p4, xdataT, ydataT, ms = 6, color = :black, label = "", inset_subplots = box,
             subplot = 2)
    # Set range of x values to plot for
    xran = 0.0:1.0:25.0
    # and then plot best fit line
    plot!(p4[2], xran, model(xran, [yintT, slopT]), label = "", color = :red, lw = wdt)
    # Set labels, fontsize, etc
    plot!(p4[2], guidefontsize = 12, legendfontsize = 12, tickfontsize = 9, grid = false)
    plot!(p4[2], xlabel = "Entropy production peaks", ylabel = "New metabolites",
          ylim = (0, 25), xlim = (0, 25))
    savefig(p4, "Output/Fig2/entp.png")
    # Now want to make a plot incorporating all four previous plots
    pt = plot(p1, p3, p2, p4, layout = (4, 1), size = (900, 1600), margin = 5mm,
              grid = false)
    savefig(pt, "Output/Fig2/figure2.eps")
    return (nothing)
end

@time figure2(1, 5, true, 61, 250, "i", 0.01, 250)
