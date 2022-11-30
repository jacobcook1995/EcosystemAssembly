# Script to construct figure 3
using Assembly
using Plots
using JLD
using StatsPlots
using StatsBase
using Statistics
using LaTeXStrings
using ColorSchemes
using Plots.PlotMeasures
using KernelDensity
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

# function to calculate simpson indices of diversity for a set of data
function simp_ind(Rl::Int64, Ru::Int64, syn::Bool, en::String, Ni::Int64, rps::Int64,
                  ipop::Float64)
    # Preallocate output data
    Sindi = zeros(rps)
    Sindf = zeros(rps)
    # Loop over repeats
    for i in 1:rps
        # Read in relevant files
        pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
        if ~isfile(pfile)
            error("run $(Nr) is missing a parameter file")
        end
        ofile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RedOutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
        if ~isfile(ofile)
            error("run $(Nr) is missing an output file")
        end
        efile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RedExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
        if ~isfile(efile)
            error("run $(Nr) is missing an extinct file")
        end
        # Load required data
        ps = load(pfile, "ps")
        inf_out = load(ofile, "inf_out")
        ded = load(efile, "ded")
        # Make containers to store initial and final types relative populations
        itype = zeros(ps.M - 1)
        ftype = zeros(ps.M - 1)
        # Loop over surviving strains for initial and final pops
        for j in 1:(ps.N)
            # Find index of main reaction
            _, ind = findmax(ps.mics[j].ϕP)
            # Find reaction
            r = ps.reacs[ps.mics[j].Reacs[ind]]
            # Add initial population to corresponding functional group
            itype[r.Rct] += ipop
            # Add relevant final populations to functional group
            ftype[r.Rct] += inf_out[j]
        end
        # Then loop over extinct strains just for initial pops
        for j in eachindex(ded)
            # Find index of main reaction
            _, ind = findmax(ded[j].ϕP)
            # Find reaction
            r = ps.reacs[ded[j].Reacs[ind]]
            # Add initial population to corresponding functional group
            itype[r.Rct] += ipop
        end
        # Normalise vectors of abundances
        itype = itype / sum(itype)
        ftype = ftype / sum(ftype)
        # Preallocate simpsons index
        λi = 0.0
        λf = 0.0
        # Loop over relative abundances to calculate simpsons index
        for j in eachindex(itype)
            λi += itype[j]^2
            λf += ftype[j]^2
        end
        Sindi[i] = 1 / λi
        Sindf[i] = 1 / λf
    end
    return (Sindi, Sindf)
end

# Function to make bar chart of the diversity loss
function divloss(Rl::Int64, Ru::Int64, syn::Bool, en::String, Ni::Int64, Nr::Int64,
                 rps::Int64)
    # Read in relevant files
    pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/ParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(Nr)Ns$(Ni).jld"
    if ~isfile(pfile)
        error("run $(Nr) is missing a parameter file")
    end
    ofile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/OutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(Nr)Ns$(Ni).jld"
    if ~isfile(ofile)
        error("run $(Nr) is missing an output file")
    end
    efile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/ExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(Nr)Ns$(Ni).jld"
    if ~isfile(efile)
        error("run $(Nr) is missing an extinct file")
    end
    # Load everything that's potentially useful
    ps = load(pfile, "ps")
    C = load(ofile, "C")
    T = load(ofile, "T")
    out = load(ofile, "out")
    ded = load(efile, "ded")
    # Remake initial list of microbes
    ms = Array{MicrobeP, 1}(undef, Ni)
    # Setup counter
    cnt = 0
    # Loop over all microbes
    for i in 1:Ni
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
    r_cat = zeros(Int64, Ni)
    # Loop over every microbe in here and categorise
    for i in 1:Ni
        # Find index of main reaction
        _, ind = findmax(ms[i].ϕP)
        # Find reaction
        r = ps.reacs[ms[i].Reacs[ind]]
        # And store index
        r_cat[i] = r.Rct
    end
    # Choose number of time points to do
    Tp = 8
    # Find and save maximum time
    Tmax = T[end]
    # Make vector of the time points
    Ts = [collect(range(0.0, stop = Tmax / 25.0, length = Tp - 1)); Tmax]
    # Container to time varying diversity and abundance data
    abT = zeros(Float64, Tp, ps.M)
    # Loop over reactions
    for i in 1:(ps.M)
        # Find all indices for this reaction
        inds = findall(x -> x == i, r_cat)
        for j in 1:Tp
            # Find first index past time point
            Tind = findfirst(x -> x >= Ts[j], T)
            # Sum abundances of survivors
            if Tind != 1
                # Interpolate abundances
                dT = ((Ts[j] - T[Tind - 1]) * sum(C[Tind, inds]) +
                      (T[Tind] - Ts[j]) * sum(C[Tind - 1, inds])) / (T[Tind] - T[Tind - 1])
                # Then save
                abT[j, i] = dT
            else
                abT[j, i] = sum(C[Tind, inds])
            end
        end
    end
    # Rescale bars to 1
    for i in 1:Tp
        abT[i, :] = abT[i, :] / (sum(abT[i, :]))
    end
    # Create xlabels
    xs = fill("", Tp)
    for i in 1:Tp
        if Ts[i] != 0.0
            # Round time
            Tr = round(Ts[i], sigdigits = 2)
            # Find power of ten
            p10 = floor(Int64, log10(Tr))
            # reduce T by this factor of 10
            rT = Tr / (10.0^p10)
            # Put together as label
            xs[i] = "$(rT)e$(p10)"
        else
            xs[i] = "0.0"
        end
    end
    # Load in colorschemes
    a = ColorSchemes.tab20.colors
    b = ColorSchemes.tab20b.colors
    # Empty array to store colours
    cls = Array{RGB{Float64}, 1}(undef, 40)
    for i in 1:20
        cls[i] = a[i]
    end
    for i in 21:40
        cls[i] = b[i - 20]
    end
    # Make my own colour scheme
    sch = ColorScheme(cls)
    # Find simpson indices
    Si, Sf = simp_ind(Rl, Ru, syn, en, Ni, rps, C[1, 1])
    # Now setup plotting
    pyplot(dpi = 200)
    p = groupedbar(abT, bar_position = :stack, label = "", palette = sch)
    plot!(p, xticks = (1:Tp, xs), ylabel = "Relative abundance", xlabel = "Time (s)")
    plot!(p, title = "Diversity with time", guidefontsize = 13, legendfontsize = 8,
          tickfontsize = 11)
    # Add annotation
    px, py = annpos([0.25; convert(Float64, Tp)], [1.0; 0.0], 0.125, 0.05)
    annotate!(px, py, text("A", 17, :black))
    # Make shape to plot beneath where the inset is going to go
    sh = Shape(5.63 .+ [0, 2.85, 2.85, 0], 0.575 .+ [0, 0, 0.135, 0.135])
    # Then plot semi-transparently
    plot!(p, sh, color = :white, lc = :white, linewidth = 0, opacity = 0.35, label = "")
    # Define box for subplot
    box = (1, bbox(0.63, 0.7, 0.31, 0.25, :bottom, :left))
    # Find kernel densities directly
    dei = kde(Si, boundary = (0, ps.M))
    def = kde(Sf, boundary = (0, ps.M))
    # Plot histogram into the subplot
    plot!(p, dei, color = :black, label = "Initial", inset_subplots = box, subplot = 2)
    plot!(p[2], def, color = :red, label = "Final", xlabel = L"^{2}D")
    # 8 and 7 before
    plot!(p[2], guidefontsize = 11, legendfontsize = 9, tickfontsize = 9, yaxis = false,
          grid = false, legend = :top)
    savefig(p, "Output/Fig3/abT.png")
    return (p)
end

function figure3(Rls::Array{Int64, 1}, Rus::Array{Int64, 1}, syns::Array{Bool, 1},
                 ens::Array{String, 1},
                 Ni::Int64, Nr::Int64, dRl::Int64, dRu::Int64, dsyn::Bool, den::String,
                 runN::Int64)
    # Check if all these vectors are the same length
    if length(Rls) != length(Rus) || length(Rls) != length(syns) ||
       length(Rls) != length(ens)
        error("length of vectors doesn't match")
    end
    println("Compiled!")
    # Count number of parameter sets
    Ns = length(Rls)
    # Container to store number of survivors
    svs = zeros(Int64, Ns, Nr)
    # Container to store metabolite diversity
    mbs = zeros(Int64, Ns, Nr)
    # Container to store dissipation
    dsp = zeros(Float64, Ns, Nr)
    # Container to store functional diversities
    fdv = zeros(Int64, Ns, Nr)
    # Loop over parameter sets
    for i in 1:Ns
        for j in 1:Nr
            # Read in relevant files
            pfile = "Data/$(Rls[i])-$(Rus[i])$(syns[i])$(Ni)$(ens[i])/RedParasReacs$(Rls[i])-$(Rus[i])Syn$(syns[i])Run$(j)Ns$(Ni).jld"
            if ~isfile(pfile)
                error("parameter set $(i) run $(j) is missing a parameter file")
            end
            ofile = "Data/$(Rls[i])-$(Rus[i])$(syns[i])$(Ni)$(ens[i])/RedOutputReacs$(Rls[i])-$(Rus[i])Syn$(syns[i])Run$(j)Ns$(Ni).jld"
            if ~isfile(ofile)
                error("parameter set $(i) run $(j) is missing an output file")
            end
            # Only want final parameter set
            ps = load(pfile, "ps")
            inf_out = load(ofile, "inf_out")
            # Save number of survivors
            svs[i, j] = ps.N
            # Loop over metabolites to find those with non-zero concentrations
            cm = 0 # Set up counter
            mm = 0 # Lowest metabolite
            for k in 2:(ps.M)
                if inf_out[ps.N + k] > 0.0
                    # Increment counter
                    cm += 1
                    # And save new minimum metabolite
                    mm = k
                end
            end
            # Save results to vector
            mbs[i, j] = cm
            # Find and save dissipation
            dsp[i, j] = dissipation(ps, ps.mics, inf_out)
            # Preallocate vector to store if function is present
            fnc = fill(false, ps.M)
            # Loop over every surviving microbe
            for k in 1:(ps.N)
                # Find index of main reaction
                _, ind = findmax(ps.mics[k].ϕP)
                # Find corresponding reactant number
                Rn = ps.reacs[ps.mics[k].Reacs[ind]].Rct
                # Update vector to show that function Rn is present
                if fnc[Rn] == false
                    fnc[Rn] = true
                end
            end
            # Count number of functions
            fdv[i, j] = count(fnc)
        end
    end
    # Setup plotting
    pyplot()
    theme(:wong2, dpi = 300)
    wongc = wong2_palette()
    # Find corresponding indices of reordered labels
    pos = zeros(Float64, Ns)
    c = Array{RGBA, 1}(undef, Ns)
    # Find positions for mean and std
    for i in 1:Ns
        p = 1.0
        if syns[i] == true
            p += 0.5
            # Set colour
            c[i] = wongc[1]
        else
            # Choose nice grey as "light" colour
            c[i] = RGB(([170, 170, 170] / 255)...)
        end
        if ens[i] == "l"
            p += 1.5
        end
        pos[i] = p
    end
    # Make latex label
    JKs = L"JK^{-1}s^{-1}"
    # Container to store mean + sd for each case
    msd = zeros(Ns)
    # Want to do the plotting here
    p1 = plot(ylabel = "Number of surviving strains", xlim = (0.5, 3.5),
              xlabel = "Energy supply")
    plot!(p1, xticks = ([1.25, 2.75], ["high", "low"]))
    # Plot means
    for i in 1:Ns
        # Calculate mean
        mn = mean(svs[i, :])
        # Calculate 99% confidence interval
        sdn = sem(svs[i, :]) * 2.576
        scatter!(p1, [pos[i]], [mn], yerror = [sdn], label = "", color = c[i], ms = 6,
                 msc = c[i])
    end
    # Add bracket for significance plot
    plot!(p1, [2.5, 3.0], [5.0, 5.0], color = :black, label = "")
    plot!(p1, [2.5, 2.5], [4.6, 5.01], color = :black, label = "")
    plot!(p1, [3.0, 3.0], [4.6, 5.01], color = :black, label = "")
    # Then add star above the bracket
    scatter!(p1, [2.75], [5.25], color = :black, shape = :star6, label = "")
    savefig(p1, "Output/Fig3/Diversity.png")
    p3 = plot(ylabel = "Entropy production rate ($(JKs))", yaxis = :log10)
    plot!(p3, xlim = (0.5, 3.5), xticks = ([1.25, 2.75], ["high", "low"]),
          xlabel = "Energy supply")
    # Plot means
    for i in 1:Ns
        # Calculate mean
        mn = mean(dsp[i, :])
        # Calculate 99% confidence interval
        sdn = sem(dsp[i, :]) * 2.576
        scatter!(p3, [pos[i]], [mn], yerror = [sdn], label = "", color = c[i], ms = 6,
                 msc = c[i])
    end
    # Add bracket for significance plot
    plot!(p3, [2.5, 3.0], [6.25e-2, 6.25e-2], color = :black, label = "")
    plot!(p3, [2.5, 2.5], [5.4e-2, 6.26e-2], color = :black, label = "")
    plot!(p3, [3.0, 3.0], [5.4e-2, 6.26e-2], color = :black, label = "")
    # Then add star above the bracket
    scatter!(p3, [2.75], [6.9e-2], color = :black, shape = :star6, label = "")
    savefig(p3, "Output/Fig3/EntropyProduction.png")
    # Want to do the plotting here
    p2 = plot(ylabel = "Survivors per substrate", xlim = (0.5, 3.5))
    plot!(p2, xticks = ([1.25, 2.75], ["high", "low"]), xlabel = "Energy supply")
    # Plot means
    for i in 1:Ns
        # Calculate mean
        mn = mean(svs[i, :] ./ mbs[i, :])
        # Calculate 99% confidence interval
        sdn = sem(svs[i, :] ./ mbs[i, :]) * 2.576
        scatter!(p2, [pos[i]], [mn], yerror = [sdn], label = "", color = c[i], ms = 6,
                 msc = c[i])
    end
    savefig(p2, "Output/Fig3/Ratio.png")
    # Want to do the plotting here
    p4 = plot(ylabel = "Number of surviving functional groups", xlim = (0.5, 3.5))
    plot!(p4, xticks = ([1.25, 2.75], ["high", "low"]), xlabel = "Energy supply",
          legend = :right)
    # Plot means
    for i in 1:Ns
        # Calculate mean
        mn = mean(fdv[i, :])
        # Calculate 99% confidence interval
        sdn = sem(fdv[i, :]) * 2.576
        # Save for annotation
        msd[i] = mn + sdn
        # Label empty
        lb = ""
        # Unless energy is high
        if ens[i] == "h"
            if syns[i] == true
                lb = "Reversible"
            else
                lb = "M–M"
            end
        end
        scatter!(p4, [pos[i]], [mn], yerror = [sdn], label = lb, color = c[i], ms = 6,
                 msc = c[i])
    end
    # Add annotation
    px, py = annpos([0.5; 3.5], msd, 0.5, -0.01)
    annotate!(px, py, text("B", 17, :black))
    # Add bracket for significance plot
    plot!(p4, [2.5, 3.0], [4.0, 4.0], color = :black, label = "")
    plot!(p4, [2.5, 2.5], [3.7, 4.01], color = :black, label = "")
    plot!(p4, [3.0, 3.0], [3.7, 4.01], color = :black, label = "")
    # Then add star above the bracket
    scatter!(p4, [2.75], [4.2], color = :black, shape = :star6, label = "")
    savefig(p4, "Output/Fig3/FuncDiv.png")
    # Combine all three plots into a single one
    pc = plot(p4, p1, p2, p3, layout = (1, 4), size = (800, 400), guidefontsize = 13,
              legendfontsize = 8, tickfontsize = 11)
    savefig(pc, "Output/Fig3/condensed.png")
    # # Run div loss function to make extra plot
    pd = divloss(dRl, dRu, dsyn, den, Ni, runN, Nr)
    # Now want to make a plot incorporating all four previous plots
    pt = plot(pd, pc, layout = grid(1, 2, widths = [0.43, 0.57]), size = (1400, 400),
              margin = 5.0mm, grid = false)
    savefig(pt, "Output/Fig3/figure3.png")
    return (nothing)
end

# Hard code parameters here
l = [1, 1, 1, 1]
u = [5, 5, 5, 5]
s = [true, true, false, false]
e = ["l", "h", "l", "h"]

@time figure3(l, u, s, e, 250, 250, 1, 5, true, "i", 89)
