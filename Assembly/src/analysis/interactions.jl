# Script to find and analyse the interaction strength and type of strains in my simulation
using Assembly
using JLD
using SymPy
using LinearAlgebra
using Statistics
using Plots
import PyPlot
using LaTeXStrings

# function to quantify the interactions
function quantify_ints()
    # Check that sufficient arguments have been provided
    if length(ARGS) < 6
        error("Insufficient inputs provided (looking for 6)")
    end
    # Preallocate the variables I want to extract from the input
    Rl = 0
    Ru = 0
    syn = true
    rps = 0
    Ni = 0
    en = ARGS[6]
    # Check that all arguments can be converted to integers
    try
        Rl = parse(Int64, ARGS[1])
        Ru = parse(Int64, ARGS[2])
        syn = parse(Bool, ARGS[3])
        rps = parse(Int64, ARGS[4])
        Ni = parse(Int64, ARGS[5])
    catch e
        error("need to provide 4 integers, a bool and a string")
    end
    # Check that simulation type is valid
    if Rl < 1
        error("lower bound has to be at least one reaction")
    end
    if Ru < Rl
        error("upper bound can't be lower than the lower bound")
    end
    # Check that number of simulations is greater than 0
    if rps < 1
        error("number of repeats can't be less than 1")
    end
    # Check that initial number of strains is sensible
    if Ni < 1
        error("number of strains cannot be less than 1")
    end
    println("Compiled!")
    flush(stdout)
    # Loop over the repeats
    for i in 1:rps
        # Assume that output files don't already exist
        outp = false
        # Just one output files to check the existence of
        outf = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/IntsReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
        # Check if all three exist
        if isfile(outf)
            outp = true
        end
        # Check if output already exists
        if outp == true
            println("Simulation $(i) already has output")
            flush(stdout)
        else
            println("Started run $(i)")
            flush(stdout)
            # Read in relevant files
            pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
            if ~isfile(pfile)
                error("run $(i) is missing a parameter file")
            end
            ofile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RedOutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
            if ~isfile(ofile)
                error("run $(i) is missing an output file")
            end
            efile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RedExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
            if ~isfile(efile)
                error("run $(i) is missing an extinct file")
            end
            # Basically just loading everything out as I'm not sure what I'll need
            ps = load(pfile, "ps")
            C = load(ofile, "C")
            T = load(ofile, "T")
            out = load(ofile, "out")
            inf_out = load(ofile, "inf_out")
            ded = load(efile, "ded")
            # For each strain record the metabolites it interacts with
            mts = Array{Array{Int64, 1}, 1}(undef, ps.N)
            # Preallocate vector of forces
            F = Array{Sym, 1}(undef, 3 * ps.N + ps.M)
            # Now find vector of forces
            F = Force(ps, F)
            # Find initial unperturbed force
            Fi = nForce(F, inf_out, ps)
            # Loop over all strains
            for j in 1:(ps.N)
                # Find substrates
                sbs = ps.reacs[ps.mics[j].Reacs] .↦ :Rct
                # Find products
                pds = ps.reacs[ps.mics[j].Reacs] .↦ :Prd
                # Combine into 1 vector (remove repeats)
                met = unique(cat(sbs, pds, dims = 1))
                # And then store
                mts[j] = met
            end
            # Preallocate vector of forces
            F = Array{Sym, 1}(undef, 3 * ps.N + ps.M)
            # Now find vector of forces
            F = Force(ps, F)
            # Preallocate perturbation
            pout = zeros(3 * ps.N + ps.M)
            # Populations, energies, and fractions remain fixed
            pout[1:(ps.N)] = inf_out[1:(ps.N)]
            pout[(ps.N + ps.M + 1):(2 * ps.N + ps.M)] = inf_out[(ps.N + ps.M + 1):(2 * ps.N + ps.M)]
            pout[(2 * ps.N + ps.M + 1):end] = inf_out[(2 * ps.N + ps.M + 1):end]
            # Preallocate array to store interaction types
            ints = zeros(Int64, ps.N, ps.N, ps.M)
            # Preallocate array to store interaction strengths
            in_str = zeros(ps.N, ps.N, ps.M)
            # Preallocate array for forces
            Fatp = zeros(ps.N, ps.M)
            # Preallocate array for net interactions
            net_in = zeros(ps.N, ps.N)
            # Set fraction to increase concentrations by
            frc = 1e-3
            # Loop over metabolites
            for j in 1:(ps.M)
                # Skip metabolites with zero concentration
                if inf_out[ps.N + j] != 0.0
                    # Make vector of +ve and -ve contributions to flux from strains
                    pvf = zeros(ps.N)
                    nvf = zeros(ps.N)
                    # Loop over strains to find values
                    for k in 1:(ps.N)
                        # Loop over strains reactions
                        for l in 1:(ps.mics[k].R)
                            # Find reaction
                            r = ps.reacs[ps.mics[k].Reacs[l]]
                            # Find amount of enzyme
                            E = Eα(inf_out[2 * ps.N + ps.M + k], ps.mics[k], l)
                            # Check if reaction has metabolite as a reactant
                            if r.Rct == j
                                q = qs(inf_out[ps.N + r.Rct], inf_out[ps.N + r.Prd], E, l,
                                       ps.mics[k], ps.T, r)
                                nvf[k] += q
                            end
                            # Or as a product
                            if r.Prd == j
                                q = qs(inf_out[ps.N + r.Rct], inf_out[ps.N + r.Prd], E, l,
                                       ps.mics[k], ps.T, r)
                                pvf[k] += q
                            end
                        end
                    end
                    # Calculate approximate perturbation size
                    dC = (pvf .- nvf) / ps.δ[j]
                    # Loop over metabolites to make perturbation
                    for k in 1:(ps.M)
                        if k != j
                            pout[ps.N + k] = inf_out[ps.N + k]
                        else
                            pout[ps.N + k] = inf_out[ps.N + k] * (1.0 + frc)
                        end
                    end
                    # Find perturbed forces
                    nF = nForce(F, pout, ps)
                    # Only interested in the forces on the ATP
                    Fatp[:, j] = nF[(ps.N + ps.M + 1):(2 * ps.N + ps.M)] .-
                                 Fi[(ps.N + ps.M + 1):(2 * ps.N + ps.M)]
                    # Preallocate positive/negative/zero vector
                    pn0 = zeros(Int64, ps.N)
                    # Preallocate supply/consumption/zero vector
                    sc0 = zeros(Int64, ps.N)
                    # Preallocate fluxes
                    fl = zeros(ps.N)
                    # Loop over the strains
                    for k in 1:(ps.N)
                        # Check whether the effect on the strain is +ve, -ve or ~0
                        # And whether strain has reaction using metabolite j
                        if Fatp[k, j] > 1e-5 && j ∈ mts[k]
                            # Positive denoted by 1
                            pn0[k] = 1
                        elseif Fatp[k, j] < -1e-5 && j ∈ mts[k]
                            # Negative denoted by -1
                            pn0[k] = -1
                        end
                        # Check net effect of strain
                        if dC[k] > 0.0
                            # Positive denoted by 1
                            sc0[k] = 1
                        elseif dC[k] < 0.0
                            # Negative denoted by -1
                            sc0[k] = -1
                        end
                    end
                    # Loop over vectors to make interaction matrix
                    for k in 1:(ps.N)
                        for l in 1:(ps.N)
                            # Check for positive cases
                            if pn0[k] == 1
                                # Crossfeeding
                                if sc0[l] == 1
                                    ints[k, l, j] = 2
                                    in_str[k, l, j] = abs(Fatp[k, j]) * abs(dC[l])
                                    # Competition
                                elseif sc0[l] == -1
                                    ints[k, l, j] = 1
                                    in_str[k, l, j] = abs(Fatp[k, j]) * abs(dC[l])
                                end
                                # And negative case
                            elseif pn0[k] == -1
                                # Competition via pollution
                                if sc0[l] == 1
                                    ints[k, l, j] = 4
                                    in_str[k, l, j] = abs(Fatp[k, j]) * abs(dC[l])
                                    # Obligate syntrophy
                                elseif sc0[l] == -1
                                    ints[k, l, j] = 3
                                    in_str[k, l, j] = abs(Fatp[k, j]) * abs(dC[l])
                                end
                            end
                        end
                        # Rescale all strengths by the population, and perturbation size
                        in_str[k, :, j] = (inf_out[k] / (ps.mics[k].ρ * ps.mics[k].MC)) *
                                          (1 / (frc * inf_out[ps.N + j])) * in_str[k, :, j]
                    end
                end
            end
            # Rescale to moles and fraction
            in_str = in_str / (NA)
            # Want to find net reactions, loop over strains
            for j in 1:(ps.N)
                for k in 1:(ps.N)
                    # Loop over every metabolite
                    for l in 1:(ps.M)
                        # Check what type of interaction (if any)
                        # effect of k on j
                        if ints[j, k, l] == 1
                            net_in[j, k] -= in_str[j, k, l]
                        elseif ints[j, k, l] == 2
                            net_in[j, k] += in_str[j, k, l]
                        elseif ints[j, k, l] == 3
                            net_in[j, k] += in_str[j, k, l]
                        elseif ints[j, k, l] == 4
                            net_in[j, k] -= in_str[j, k, l]
                        end
                    end
                end
            end
            # Output all interaction data
            jldopen("Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/IntsReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld",
                    "w") do file
                # Save ATP forces and fraction used to generate them
                write(file, "Fatp", Fatp)
                write(file, "frc", frc)
                # Save interactions and interaction strengths
                write(file, "ints", ints)
                write(file, "in_str", in_str)
                # Save net interactions
                write(file, "net_in", net_in)
            end
        end
    end
    return (nothing)
end

# Function to analyse interaction data
function analyse_ints()
    # Check that sufficient arguments have been provided
    if length(ARGS) < 4
        error("Insufficient inputs provided (looking for 4)")
    end
    # Preallocate the variables I want to extract from the input
    Rl = 0
    Ru = 0
    syn = true
    rps = 0
    # Check that all arguments can be converted to integers
    try
        Rl = parse(Int64, ARGS[1])
        Ru = parse(Int64, ARGS[2])
        syn = parse(Bool, ARGS[3])
        rps = parse(Int64, ARGS[4])
    catch e
        error("need to provide 3 integers and a bool")
    end
    # Check that simulation type is valid
    if Rl < 1
        error("lower bound has to be at least one reaction")
    end
    if Ru < Rl
        error("upper bound can't be lower than the lower bound")
    end
    # Check that number of simulations is greater than 0
    if rps < 1
        error("number of repeats can't be less than 1")
    end
    println("Compiled!")
    # Set time to run perturbed dynamics for
    Tmax = 1e5
    # Set fraction to perturb (reduce) population by
    pfrc = 0.5
    # Setup counters
    s1 = 0 # Self interaction
    st = 0 # Correct peak for self interaction
    n1 = 0 # no first order
    n1t = 0 # correct peak (no first order)
    n2 = 0 # no second order
    n2t = 0 # correct peak (no second order)
    cc = 0 # correctly predicted consistent peak
    ci = 0 # incorrectly predicted consistent peak
    ic = 0 # correctly predicted inconsistent peak
    ii = 0 # incorrectly predicted inconsistent peak
    Tt = 2500.0 # threshold T value
    # Set up plotting
    for i in 1:rps
        # Read in relevant files
        pfile = "Data/$(Rl)-$(Ru)$(syn)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        ofile = "Data/$(Rl)-$(Ru)$(syn)/RedOutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        if ~isfile(ofile)
            error("run $(i) is missing an output file")
        end
        efile = "Data/$(Rl)-$(Ru)$(syn)/RedExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        if ~isfile(efile)
            error("run $(i) is missing an extinct file")
        end
        ifile = "Data/$(Rl)-$(Ru)$(syn)/IntsReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        if ~isfile(ifile)
            error("run $(i) is missing an interaction file")
        end
        # Basically just loading everything out as I'm not sure what I'll need
        ps = load(pfile, "ps")
        C = load(ofile, "C")
        T = load(ofile, "T")
        out = load(ofile, "out")
        inf_out = load(ofile, "inf_out")
        ded = load(efile, "ded")
        Fatp = load(ifile, "Fatp")
        frc = load(ifile, "frc")
        ints = load(ifile, "ints")
        in_str = load(ifile, "in_str")
        net_in = load(ifile, "net_in")
        # Preallocate self-interaction and interaction matrices
        R = zeros(ps.N, ps.N)
        G = zeros(ps.N, ps.N)
        # Populate these values
        for j in 1:(ps.N)
            for k in 1:(ps.N)
                if j == k
                    R[j, k] = net_in[j, k]
                else
                    G[j, k] = net_in[j, k]
                end
            end
        end
        # Make squared matrix
        G2 = G * G
        G3 = G2 * G
        G4 = G3 * G
        G5 = G4 * G
        # concentrations and internal cell parameters are not perturbed
        conc = inf_out[(ps.N + 1):(ps.N + ps.M)]
        as = inf_out[(ps.N + ps.M + 1):(2 * ps.N + ps.M)]
        ϕs = inf_out[(2 * ps.N + ps.M + 1):end]
        # Only populations are perturbed
        pop = zeros(ps.N)
        # Preallocate interaction matrix
        pet_in = zeros(ps.N, ps.N)
        # Preallocate array to store peak signs
        psg = zeros(Int64, ps.N, ps.N)
        # Perturb each strain in term
        for j in 1:(ps.N)
            # Loop over strains to make perturbation
            for k in 1:(ps.N)
                if k != j
                    pop[k] = inf_out[k]
                else
                    pop[k] = inf_out[k] * (1.0 - pfrc)
                end
            end
            # Then run the simulation
            Cp, Tp = full_simulate(ps, Tmax, pop, conc, as, ϕs)
            # Preallocate gradients and peak indices
            grds = zeros(length(Tp) - 1, ps.N)
            pinds = zeros(Int64, ps.N)
            # Find gradients for each strain
            for k in 1:(ps.N)
                for l in 2:length(Tp)
                    grds[l - 1, k] = (Cp[l, k] - Cp[l - 1, k]) / (Tp[l] - Tp[l - 1])
                end
            end
            # Find index of peak for each strain
            for k in 1:(ps.N)
                # Find largest deviation from initial value
                _, pinds[k] = findmax(abs.(Cp[:, k] .- Cp[1, k]))
                # Then calculate the sign of this peak
                psg[k, j] = sign(Cp[pinds[k], k] - Cp[1, k])
            end
            # Set range to pick value from
            rn = 1:3
            for k in 1:(ps.N)
                # Find first non-zero element
                x1 = findfirst(x -> x != 0.0, grds[:, k])
                # and use to update range
                rn2 = rn .+ (x1 - 1)
                # Find absolute maximum gradient from the range
                _, gin = findmax(abs.(grds[rn2, k]))
                # Add first element of the range to returned index
                gind = gin + (rn2[1] - 1)
                # Save gradients as "effective" interaction strengths
                pet_in[k, j] = grds[gind, k]
            end
        end
        # Rescale interactions to see if they match
        if ps.N != 0
            pin = -pet_in / (maximum(abs.(pet_in)))
            nin = net_in / (maximum(abs.(net_in)))
        end
        # Loop over all entries
        for j in 1:(ps.N)
            for k in 1:(ps.N)
                # Check if signs match
                if sign(pin[j, k]) != sign(nin[j, k]) && abs(nin[j, k]) > 5e-4
                    println("Run $(i)")
                    println("Predicted initial effect of strain $(k) on strain $(j) incorrect")
                end
                # Find time spans for each interaction
                if G2[j, k] != 0
                    T2 = abs((R[j, k] + G[j, k]) / G2[j, k])
                else
                    T2 = Inf
                end
                if G3[j, k] != 0
                    T3 = (abs((R[j, k] + G[j, k]) / G3[j, k]))^(1 / 2)
                else
                    T3 = Inf
                end
                if G4[j, k] != 0
                    T4 = (abs((R[j, k] + G[j, k]) / G4[j, k]))^(1 / 3)
                else
                    T4 = Inf
                end
                if G5[j, k] != 0
                    T5 = (abs((R[j, k] + G[j, k]) / G5[j, k]))^(1 / 4)
                else
                    T5 = Inf
                end
                # Use T values to select dominant higher order interaction
                if T2 < T3 && T2 < T4 && T2 < T5
                    GA = G2[j, k]
                    TA = T2
                elseif T3 < T4 && T3 < T5
                    GA = G3[j, k]
                    TA = T3
                elseif T4 < T5
                    GA = G4[j, k]
                    TA = T4
                else
                    GA = G5[j, k]
                    TA = T5
                end
                # Check if it's a recovery from perturbation
                if j == k
                    s1 += 1
                    # Check peak is of the expected type
                    if sign(psg[j, k]) == 1
                        st += 1
                    end
                    # Check if no interaction
                elseif R[j, k] + G[j, k] == 0.0
                    n1 += 1
                    # Check that sign of peak matches higher order interaction
                    if sign(GA) == -sign(psg[j, k])
                        n1t += 1
                    end
                    # Check if no mediated interaction
                elseif GA == 0.0
                    n2 += 1
                    # Check that peak matches initial gradient
                    if sign(nin[j, k]) == -sign(psg[j, k])
                        n2t += 1
                    end
                    # Check if higher order interaction is to slow or if interaction signs match
                elseif TA > Tt || sign(R[j, k] + G[j, k]) == sign(GA)
                    if sign(nin[j, k]) != -sign(psg[j, k])
                        ci += 1
                    else
                        cc += 1
                    end
                else
                    if sign(nin[j, k]) != sign(psg[j, k])
                        ii += 1
                    else
                        ic += 1
                    end
                end
            end
        end
    end
    println("$(s1) recoveries from perturbations in total")
    println("$(st) actual recoveries observed")
    println("$(n1+n2+cc+ci+ic+ii) peaks in total")
    println("$(n1) peaks with no 1st order interaction")
    println("Of which $(n1t) are consistent")
    println("$(n2) peaks with no 2nd order interaction")
    println("Of which $(n2t) are consistent")
    println("$(cc+ci) consistent peaks predicted")
    println("$(cc) consistent peaks found")
    println("$(ic+ii) inconsistent peaks predicted")
    println("$(ic) inconsistent peaks found")
    println("$(n1t+n2t+cc+ic) out of $(n1+n2+cc+ci+ic+ii) guessed correctly")
    return (nothing)
end

# Function to analyse interaction data
function plot_ptrbs()
    # Check that sufficient arguments have been provided
    if length(ARGS) < 4
        error("Insufficient inputs provided (looking for 4)")
    end
    # Preallocate the variables I want to extract from the input
    Rl = 0
    Ru = 0
    syn = true
    rpt = 0
    # Check that all arguments can be converted to integers
    try
        Rl = parse(Int64, ARGS[1])
        Ru = parse(Int64, ARGS[2])
        syn = parse(Bool, ARGS[3])
        rpt = parse(Int64, ARGS[4])
    catch e
        error("need to provide 3 integers and a bool")
    end
    # Check that simulation type is valid
    if Rl < 1
        error("lower bound has to be at least one reaction")
    end
    if Ru < Rl
        error("upper bound can't be lower than the lower bound")
    end
    # Check that number of simulations is greater than 0
    if rpt < 1
        error("number of repeats can't be less than 1")
    end
    println("Compiled!")
    # Set time to run perturbed dynamics for
    Tmax = 1e5
    # Set fraction to perturb (reduce) population by
    pfrc = 0.5
    pyplot()
    theme(:wong2, dpi = 200)
    wongc = get_color_palette(wong_palette, 57)
    # Set up plotting
    for i in 1:rpt
        # Read in relevant files
        pfile = "Data/$(Rl)-$(Ru)$(syn)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        ofile = "Data/$(Rl)-$(Ru)$(syn)/RedOutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        if ~isfile(ofile)
            error("run $(i) is missing an output file")
        end
        efile = "Data/$(Rl)-$(Ru)$(syn)/RedExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        if ~isfile(efile)
            error("run $(i) is missing an extinct file")
        end
        ifile = "Data/$(Rl)-$(Ru)$(syn)/IntsReacs$(Rl)-$(Ru)Syn$(syn)Run$(i).jld"
        if ~isfile(ifile)
            error("run $(i) is missing an interaction file")
        end
        # Basically just loading everything out as I'm not sure what I'll need
        ps = load(pfile, "ps")
        C = load(ofile, "C")
        T = load(ofile, "T")
        out = load(ofile, "out")
        inf_out = load(ofile, "inf_out")
        ded = load(efile, "ded")
        Fatp = load(ifile, "Fatp")
        frc = load(ifile, "frc")
        ints = load(ifile, "ints")
        in_str = load(ifile, "in_str")
        net_in = load(ifile, "net_in")
        # concentrations and internal cell parameters are not perturbed
        conc = inf_out[(ps.N + 1):(ps.N + ps.M)]
        as = inf_out[(ps.N + ps.M + 1):(2 * ps.N + ps.M)]
        ϕs = inf_out[(2 * ps.N + ps.M + 1):end]
        # Only populations are perturbed
        pop = zeros(ps.N)
        # Perturb each strain in term
        for j in 1:(ps.N)
            # Loop over strains to make perturbation
            for k in 1:(ps.N)
                if k != j
                    pop[k] = inf_out[k]
                else
                    pop[k] = inf_out[k] * (1.0 - pfrc)
                end
            end
            # Then run the simulation
            Cp, Tp = full_simulate(ps, Tmax, pop, conc, as, ϕs)
            # Preallocate gradients and peak indices
            grds = zeros(length(Tp) - 1, ps.N)
            pinds = zeros(Int64, ps.N)
            # Find gradients for each strain
            for k in 1:(ps.N)
                for l in 2:length(Tp)
                    grds[l - 1, k] = (Cp[l, k] - Cp[l - 1, k]) / (Tp[l] - Tp[l - 1])
                end
            end
            # Preallocate plot
            plot(title = "Perturbed strain $(j)", yaxis = :log10, xlabel = "Time",
                 ylabel = "Log population")
            # Loop over and plot all strains
            for k in 1:(ps.N)
                # Find and eliminate zeros so that they can be plotted on a log plot
                inds = (Cp[:, k] .> 0)
                plot!(Tp[inds], Cp[inds, k], color = wongc[k], label = "Strain $(k)")
            end
            # Save figure
            savefig("Output/perturbpops$(j).png")
        end
    end
    return (nothing)
end

# function to plot interactions types
function ints_plot()
    # Check that sufficient arguments have been provided
    if length(ARGS) < 6
        error("Insufficient inputs provided (looking for 6)")
    end
    # Preallocate the variables I want to extract from the input
    Rl = 0
    Ru = 0
    syn = true
    rps = 0
    Ni = 0
    en = ARGS[6]
    # Check that all arguments can be converted to integers
    try
        Rl = parse(Int64, ARGS[1])
        Ru = parse(Int64, ARGS[2])
        syn = parse(Bool, ARGS[3])
        rps = parse(Int64, ARGS[4])
        Ni = parse(Int64, ARGS[5])
    catch e
        error("need to provide 4 integers and a bool")
    end
    # Check that simulation type is valid
    if Rl < 1
        error("lower bound has to be at least one reaction")
    end
    if Ru < Rl
        error("upper bound can't be lower than the lower bound")
    end
    # Check that number of simulations is greater than 0
    if rps < 1
        error("number of repeats can't be less than 1")
    end
    if Ni < 1
        error("initial number of strains can't be less than 1")
    end
    println("Compiled!")
    # Preallocate memory to store number of interactions
    ins1 = zeros(rps)
    ins2 = zeros(rps)
    ins3 = zeros(rps)
    ins4 = zeros(rps)
    # Preallocate memory to store mean interaction strengths
    mean1 = fill(NaN, rps)
    mean2 = fill(NaN, rps)
    mean3 = fill(NaN, rps)
    mean4 = fill(NaN, rps)
    # Preallocate vectors to store all interaction strengths
    sts1 = []
    sts2 = []
    sts3 = []
    sts4 = []
    # Loop over parameter sets
    for i in 1:rps
        # Read in relevant files
        pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        ofile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RedOutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
        if ~isfile(ofile)
            error("run $(i) is missing an output file")
        end
        efile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RedExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
        if ~isfile(efile)
            error("run $(i) is missing an extinct file")
        end
        ifile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/IntsReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
        if ~isfile(ifile)
            error("run $(i) is missing an interaction file")
        end
        # Basically just loading everything out as I'm not sure what I'll need
        ps = load(pfile, "ps")
        C = load(ofile, "C")
        T = load(ofile, "T")
        out = load(ofile, "out")
        inf_out = load(ofile, "inf_out")
        ded = load(efile, "ded")
        Fatp = load(ifile, "Fatp")
        frc = load(ifile, "frc")
        ints = load(ifile, "ints")
        in_str = load(ifile, "in_str")
        net_in = load(ifile, "net_in")
        # Store number of each interaction type
        ins1[i] = count(x -> x == 1, ints)
        ins2[i] = count(x -> x == 2, ints)
        ins3[i] = count(x -> x == 3, ints)
        ins4[i] = count(x -> x == 4, ints)
        # Find corresponding positions
        pos1 = findall(x -> x == 1, ints)
        pos2 = findall(x -> x == 2, ints)
        pos3 = findall(x -> x == 3, ints)
        pos4 = findall(x -> x == 4, ints)
        # Use to find vector of interaction strengths
        str1 = in_str[pos1]
        str2 = in_str[pos2]
        str3 = in_str[pos3]
        str4 = in_str[pos4]
        # Check if there are any interactions of that type
        if length(str1) >= 1
            # Find mean strength of each interaction
            mean1[i] = mean(str1)
            sts1 = cat(sts1, str1, dims = 1)
        end
        # Do for each interaction type
        if length(str2) >= 1
            mean2[i] = mean(str2)
            sts2 = cat(sts2, str2, dims = 1)
        end
        if length(str3) >= 1
            mean3[i] = mean(str3)
            sts3 = cat(sts3, str3, dims = 1)
        end
        if length(str4) >= 1
            mean4[i] = mean(str4)
            sts4 = cat(sts4, str4, dims = 1)
        end
    end
    # Setup plotting
    pyplot()
    theme(:wong2, dpi = 200)
    wongc = get_color_palette(wong_palette, 57)
    # Make plot title
    tl = ""
    if syn == true
        tl = "$(Rl)-$(Ru) reactions per strain"
    else
        tl = "$(Rl)-$(Ru) reactions per strain (no syntrophy)"
    end
    # Now plot both interactions types and strengths
    plot(title = tl, xlabel = "Number of interactions", ylabel = "Number of ecosystems")
    histogram!(ins1, fillalpha = 0.75, label = "Competition")
    histogram!(ins2, fillalpha = 0.75, label = "Facilitation")
    histogram!(ins3, fillalpha = 0.75, label = "Syntrophy")
    histogram!(ins4, fillalpha = 0.75, label = "Pollution")
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/IntType$(Rl)-$(Ru)$(syn)$(Ni)$(en).png")
    plot(title = tl, xlabel = "Strength of interactions", ylabel = "Number of ecosystems")
    histogram!(log10.(mean1), fillalpha = 0.75, label = "Competition")
    histogram!(log10.(mean2), fillalpha = 0.75, label = "Facilitation")
    histogram!(log10.(mean3), fillalpha = 0.75, label = "Syntrophy")
    histogram!(log10.(mean4), fillalpha = 0.75, label = "Pollution")
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/IntStrength$(Rl)-$(Ru)$(syn)$(Ni)$(en).png")
    # Make range of ticks to label
    rgn = collect(-14:2:-2)
    ergn = fill("", length(rgn))
    # Find corresponding exponentials
    for i in 1:length(rgn)
        ergn[i] = L"10^{%$(rgn[i])}"
    end
    plot(title = tl, xlabel = "Interaction strength", ylabel = "Number of interactions")
    plot!(xticks = (rgn, ergn))
    histogram!(log10.(sts1), fillalpha = 0.75, label = "Competition")
    histogram!(log10.(sts2), fillalpha = 0.75, label = "Facilitation")
    histogram!(log10.(sts3), fillalpha = 0.75, label = "Syntrophy")
    histogram!(log10.(sts4), fillalpha = 0.75, label = "Pollution")
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/AllIntStrength$(Rl)-$(Ru)$(syn)$(Ni)$(en).png")
    return (nothing)
end

# function to plot interactions types
function ints_scat()
    # Check that sufficient arguments have been provided
    if length(ARGS) < 6
        error("Insufficient inputs provided (looking for 6)")
    end
    # Preallocate the variables I want to extract from the input
    Rl = 0
    Ru = 0
    syn = true
    rps = 0
    Ni = 0
    en = ARGS[6]
    # Check that all arguments can be converted to integers
    try
        Rl = parse(Int64, ARGS[1])
        Ru = parse(Int64, ARGS[2])
        syn = parse(Bool, ARGS[3])
        rps = parse(Int64, ARGS[4])
        Ni = parse(Int64, ARGS[5])
    catch e
        error("need to provide 4 integers and a bool")
    end
    # Check that simulation type is valid
    if Rl < 1
        error("lower bound has to be at least one reaction")
    end
    if Ru < Rl
        error("upper bound can't be lower than the lower bound")
    end
    # Check that number of simulations is greater than 0
    if rps < 1
        error("number of repeats can't be less than 1")
    end
    if Ni < 1
        error("initial number of strains can't be less than 1")
    end
    println("Compiled!")
    # Preallocate memory to store number of interactions
    ins3 = zeros(rps)
    # preallocate memory for syntrophic interaction strengths
    i3str = fill(NaN, rps)
    # Preallocate memory to store mean interaction strengths
    mean1 = zeros(rps)
    mean2 = zeros(rps)
    mean3 = zeros(rps)
    mean4 = zeros(rps)
    # Preallocate memory to store number of survivors
    dvs = zeros(rps)
    # Preallocate memory to store relative number of syntrophy interactions
    rsyn = zeros(rps)
    # Preallocate memory to store relative number of thermodynamic interactions
    rthm = zeros(rps)
    # Preallocate memory to store relative number of cooperative interactions
    rcp = zeros(rps)
    # Loop over parameter sets
    for i in 1:rps
        # Read in relevant files
        pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        ofile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RedOutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
        if ~isfile(ofile)
            error("run $(i) is missing an output file")
        end
        efile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RedExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
        if ~isfile(efile)
            error("run $(i) is missing an extinct file")
        end
        ifile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/IntsReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
        if ~isfile(ifile)
            error("run $(i) is missing an interaction file")
        end
        # Basically just loading everything out as I'm not sure what I'll need
        ps = load(pfile, "ps")
        inf_out = load(ofile, "inf_out")
        Fatp = load(ifile, "Fatp")
        frc = load(ifile, "frc")
        ints = load(ifile, "ints")
        in_str = load(ifile, "in_str")
        net_in = load(ifile, "net_in")
        # Save number of survivors
        dvs[i] = ps.N
        # Count number of each interaction type
        i1 = count(x -> x == 1, ints)
        i2 = count(x -> x == 2, ints)
        i3 = count(x -> x == 3, ints)
        i4 = count(x -> x == 4, ints)
        # Save number of syntrophic interactions
        ins3[i] = i3
        # Save relative number of interactions
        rsyn[i] = i3 / (i1 + i2 + i3 + i4)
        rthm[i] = (i3 + i4) / (i1 + i2 + i3 + i4)
        rcp[i] = (i3 + i2) / (i1 + i2 + i3 + i4)
        # Find corresponding positions
        pos1 = findall(x -> x == 1, ints)
        pos2 = findall(x -> x == 2, ints)
        pos3 = findall(x -> x == 3, ints)
        pos4 = findall(x -> x == 4, ints)
        # Use to find vector of interaction strengths
        str1 = in_str[pos1]
        str2 = in_str[pos2]
        str3 = in_str[pos3]
        str4 = in_str[pos4]
        # Check if there are any interactions of that type
        if length(str1) >= 1
            mean1[i] = mean(str1)
        end
        if length(str2) >= 1
            mean2[i] = mean(str2)
        end
        if length(str3) >= 1
            mean3[i] = mean(str3)
            i3str[i] = mean3[i]
        end
        if length(str4) >= 1
            mean4[i] = mean(str4)
        end
    end
    # Setup plotting
    pyplot()
    theme(:wong2, dpi = 200)
    wongc = get_color_palette(wong_palette, 57)
    # Plot diversity against number of syntrophic interactions
    scatter(ins3, dvs, label = "", xlabel = "Number of syntrophic interactions",
            ylabel = "Surviving strains")
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/NumbInts$(Rl)-$(Ru)$(syn)$(Ni).png")
    # or maybe against mean interaction strength
    plot(xscale = :log10)
    scatter!(i3str, dvs, label = "", xlabel = "Mean syntrophic interaction strength",
             ylabel = "Surviving strains")
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/StrInts$(Rl)-$(Ru)$(syn)$(Ni).png")
    # Plot relative number of syntrophy interactions against diversity
    scatter(rsyn, dvs, label = "", xlabel = "Relative number of syntrophic interactions",
            ylabel = "Surviving strains")
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RelNumbInts$(Rl)-$(Ru)$(syn)$(Ni).png")
    # Plot relative number of thermodynamic interactions against diversity
    scatter(rthm, dvs, label = "", xlabel = "Relative number of thermodynamic interactions",
            ylabel = "Surviving strains")
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RelThermInts$(Rl)-$(Ru)$(syn)$(Ni).png")
    # Plot relative number of cooperative interactions against diversity
    scatter(rcp, dvs, label = "", xlabel = "Relative number of cooperative interactions",
            ylabel = "Surviving strains")
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RelCoopInts$(Rl)-$(Ru)$(syn)$(Ni).png")
    # Plot relative strength of syntrophic interactions
    scatter(mean3 ./ (mean1 .+ mean2 .+ mean3 .+ mean4), dvs, label = "")
    plot!(xlabel = "Relative strength of syntrophic interactions",
          ylabel = "Surviving strains")
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RelSynStr$(Rl)-$(Ru)$(syn)$(Ni).png")
    # Plot relative strength of syntrophic interactions
    scatter((mean3 .+ mean4) ./ (mean1 .+ mean2 .+ mean3 .+ mean4), dvs, label = "")
    plot!(xlabel = "Relative strength of thermodynamic interactions",
          ylabel = "Surviving strains")
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RelThermStr$(Rl)-$(Ru)$(syn)$(Ni).png")
    # Plot relative strength of cooperative interactions
    scatter((mean2 .+ mean4) ./ (mean1 .+ mean2 .+ mean3 .+ mean4), dvs, label = "")
    plot!(xlabel = "Relative strength of cooperative interactions",
          ylabel = "Surviving strains")
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RelCoopStr$(Rl)-$(Ru)$(syn)$(Ni).png")
    return (nothing)
end

@time ints_plot()
