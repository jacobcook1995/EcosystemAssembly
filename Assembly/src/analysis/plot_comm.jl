# This script exists to plot the output of form_comm.jl
using Assembly
using Plots
using LaTeXStrings
using JLD
import PyPlot

# function to calculate the dissipation for an assembled ecosystem
function dissipation(ps::FullParameters,ms::Array{MicrobeP,1},out::Array{Float64,1})
    # Define number of strains
    N = length(ms)
    # check that parameter set is sensible given the output
    if length(out) != ps.M + 3*N
        error("parameter set doesn't match output")
    end
    # Set dissipation to zero
    dsp = 0
    # Loop over number of strains
    for i = 1:N
        # Isolate this strain
        mic = ms[i]
        # Loop over reactions of this strain
        for j = 1:mic.R
            # Find appropriate reaction
            r = ps.reacs[mic.Reacs[j]]
            # If there's no product Gibbs free energy becomes infinite. Justified to ignore
            # this as if product hasn't built up reaction can't be happening to a sigificant degree
            if out[N+r.Prd] != 0.0
                # Find amount of energy that this reaction dissipates
                Fd = -(r.ΔG0 + Rgas*ps.T*log(out[N+r.Prd]/out[N+r.Rct]) + mic.η[j]*ΔGATP)
                # Find amount of enzyme E
                E = Eα(out[2*N+ps.M+i],mic,j)
                # Then find the rate that this reaction proceeds at
                q = qs(out[N+r.Rct],out[N+r.Prd],E,j,mic,ps.T,r)
                # Check if reaction actually occurs
                if q != 0.0
                    dsp += q*Fd*out[i]
                end
            end
        end
    end
    # Convert from molecule units to moles
    dsp /= NA
    return(dsp)
end

# A function to plot the stability time graphs
function stab_plots()
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
    # Now move onto plotting
    pyplot(dpi=200)
    # Pick random number to select parameter set
    r = rand(1:rps)
    # Read in relevant files
    pfile = "Data/Type$(R)/ParasType$(R)Run$(r).jld"
    if ~isfile(pfile)
        error("run $(r) is missing a parameter file")
    end
    ofile = "Data/Type$(R)/OutputType$(R)Run$(r).jld"
    if ~isfile(ofile)
        error("run $(r) is missing an output file")
    end
    efile = "Data/Type$(R)/ExtinctType$(R)Run$(r).jld"
    if ~isfile(efile)
        error("run $(r) is missing an extinct file")
    end
    ps = load(pfile,"ps")
    C = load(ofile,"C")
    T = load(ofile,"T")
    out = load(ofile,"out")
    ded = load(efile,"ded")
    # Setup population plot
    p1 = plot(xlabel="Time",ylabel="Population",yaxis=:log10,title="$(R) reactions")
    # Find number of ignored microbes
    N = ps.N + length(ded)
    for i = 1:N
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0)
        plot!(p1,T[inds],C[inds,i],label="")
    end
    savefig(p1,"Output/SeriesType$(R).png")
    # Now plot total population
    plot(xlabel="Time",ylabel="Total population",yaxis=:log10,title="$(R) reactions")
    plot!(T,sum(C[:,1:N],dims=2),label="")
    savefig("Output/SeriesType$(R)Total.png")
    # Preallocate vector of total population growth rates
    λs = zeros(length(T))
    # Now find growth rates for total population
    for i = 2:length(T)
        λs[i] = abs((sum(C[i,1:N]) - sum(C[i-1,1:N]))/(T[i] - T[i-1]))
        # Remove zero entries
        if λs[i] == 0.0
            λs[i] = 1e-10
        end
    end
    # Now plot growth rate of total populations
    plot(xlabel="Time",ylabel="Population change (cells/s)",yaxis=:log10,title="$(R) reactions")
    plot!(T[2:end],λs[2:end],label="")
    savefig("Output/SeriesType$(R)Pcnt.png")
    # First time frame 0 to 1e5
    p1 = plot(xlabel="Time",ylabel="Population",yaxis=:log10,title="$(R) reactions")
    p2 = plot(xlabel="Time",ylabel="Total Population",yaxis=:log10,title="$(R) reactions")
    p3 = plot(xlabel="Time",ylabel="Population change (cells/s)",yaxis=:log10,title="$(R) reactions")
    # Find number of ignored microbes
    N = ps.N + length(ded)
    for i = 1:N
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0) .& (T .< 5e5)
        plot!(p1,T[inds],C[inds,i],label="")
    end
    # Find indices of points less than particular time
    inds = (T .< 5e5)
    plot!(p2,T[inds],sum(C[inds,:],dims=2),label="")
    # Change indicesso that the first point is ignored
    inds2 = inds
    inds2[1] = false
    plot!(p3,T[inds],λs[inds],label="")
    savefig(p1,"Output/SeriesType$(R)F1.png")
    savefig(p2,"Output/SeriesType$(R)TotalF1.png")
    savefig(p3,"Output/SeriesType$(R)PcntF1.png")
    # Second time frame 0 to 1e7
    p1 = plot(xlabel="Time",ylabel="Population",yaxis=:log10,title="$(R) reactions")
    p2 = plot(xlabel="Time",ylabel="Total Population",yaxis=:log10,title="$(R) reactions")
    p3 = plot(xlabel="Time",ylabel="Population change (cells/s)",yaxis=:log10,title="$(R) reactions")
    # Find number of ignored microbes
    N = ps.N + length(ded)
    for i = 1:N
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0) .& (T .< 5e6)
        plot!(p1,T[inds],C[inds,i],label="")
    end
    # Find indices of points less than particular time
    inds = (T .< 5e6)
    plot!(p2,T[inds],sum(C[inds,:],dims=2),label="")
    # Change indices that the first point is ignored
    inds2 = inds
    inds2[1] = false
    plot!(p3,T[inds],λs[inds],label="")
    savefig(p1,"Output/SeriesType$(R)F2.png")
    savefig(p2,"Output/SeriesType$(R)TotalF2.png")
    savefig(p3,"Output/SeriesType$(R)PcntF2.png")
    # Third time frame 5e7 to maximum
    p1 = plot(xlabel="Time",ylabel="Population",yaxis=:log10,title="$(R) reactions")
    p2 = plot(xlabel="Time",ylabel="Total Population",yaxis=:log10,title="$(R) reactions")
    p3 = plot(xlabel="Time",ylabel="Population change (cells/s)",yaxis=:log10,title="$(R) reactions")
    # Find number of ignored microbes
    N = ps.N + length(ded)
    for i = 1:N
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0) .& (T .> 5e7)
        plot!(p1,T[inds],C[inds,i],label="")
    end
    # Find indices of points greater than particular time
    inds = (T .> 5e7)
    plot!(p2,T[inds],sum(C[inds,:],dims=2),label="")
    # Change indices that the first point is ignored
    inds2 = inds
    inds2[1] = false
    plot!(p3,T[inds],λs[inds],label="")
    savefig(p1,"Output/SeriesType$(R)F3.png")
    savefig(p2,"Output/SeriesType$(R)TotalF3.png")
    savefig(p3,"Output/SeriesType$(R)PcntF3.png")
    # Preallocate vector of extinctions
    ext = BitArray(undef,ps.N)
    # Loop over every microbe
    for i = 1:ps.N
        # Find microbes index in full dynamics
        ind = findfirst(x->x==out[i],C[end,:])
        # Calculate rate change in final two steps of the dynamics
        dC = ((C[end,ind] - C[end-1,ind])/(T[end] - T[end-1]))/(C[end,ind])
        # If strain decreasing at greater than threshold rate then assume extinct
        if dC <= -5e-9
            ext[i] = 1
        else
            ext[i] = 0
        end
    end
    # Now need to find reduced parameter set
    ps = extinction(ps,ext)
    # Needs to run for a good while to confirm stability
    Tmax = 1e7
    # Preallocate vectors to store the intial conditions
    pop = zeros(ps.N)
    as = zeros(ps.N)
    ϕs = zeros(ps.N)
    # Concentration data can be imediatly filled in
    conc = out[length(ext)+1:length(ext)+ps.M]
    # counter function
    k = 0
    # Now fill in the data
    for i = 1:length(ext)
        if ext[i] == 0
            pop[i-k] = out[i]
            as[i-k] = out[length(ext)+ps.M+i]
            ϕs[i-k] = out[2*length(ext)+ps.M+i]
        else
            k += 1
        end
    end
    # Run the simulation
    C, T = full_simulate(ps,Tmax,pop,conc,as,ϕs)
    # Now do standard population plot
    p1 = plot(xlabel="Time",ylabel="Population",yaxis=:log10,title="Stablity check $(R) reactions")
    for i = 1:ps.N
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0)
        plot!(p1,T[inds],C[inds,i],label="")
    end
    savefig("Output/StabCheckType$(R).png")
    return(nothing)
end

# New function to plot histograms over time
function hist_time()
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
    # Choose time steps
    Ts = [1e4,1e5,1e6,1e7,1e8]
    # Preallocate memory to store the number of survivors in
    srv = zeros(rps,length(Ts))
    # Memory to store number of expected survivors (at infinity)
    svinf = zeros(rps,length(Ts))
    # Preallocate vector to store percentage of free energy dissipated (under standard conditions)
    pd = fill(Float64[],length(Ts))
    # Similar vector to store populations
    pps = fill(Float64[],length(Ts))
    # Preallocate vector of dissipation rates
    dsp = zeros(rps,length(Ts))
    # And the same for dissipation rate per unit biomass
    b_dsp = zeros(rps,length(Ts))
    # Preallocate catagories
    inact = zeros(rps,length(Ts))
    ob = zeros(rps,length(Ts))
    nob = zeros(rps,length(Ts))
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
        if ~isfile(efile)
            error("run $(i) is missing an extinct file")
        end
        # Then load in data
        ps = load(pfile,"ps")
        out = load(ofile,"out")
        T = load(ofile,"T")
        C = load(ofile,"C")
        ded = load(efile,"ded")
        # Find number of intial strains
        N = ps.N + length(ded)
        # Preallocate vector of microbes
        ms = Vector{MicrobeP}(undef,N)
        # Find indices of surviving microbes
        inds = indexin(out[1:ps.N],C[end,1:N])
        # Set up loop over all microbes
        k = 0
        for j = 1:length(ms)
            # Add surviving strains in
            if j ∈ inds
                ms[j] = ps.mics[j-k]
            # and extinct strains
            else
                k += 1
                ms[j] = ded[k]
            end
        end
        # Loop over time steps
        for j = 1:length(Ts)
            # Find index of first point that is a greater time
            indT = findfirst(x->x>=Ts[j],T)
            # Find distance between this point and the prior point
            dT = T[indT] - T[indT-1]
            # Use this to find weight this time point should be given
            w = 1 - (T[indT] - Ts[j])/dT
            # Make new vectors for populations and metabolites
            pops = w*C[indT,1:N] .+ (1-w)*C[indT-1,1:N]
            concs = w*C[indT,N+1:N+ps.M] .+ (1-w)*C[indT-1,N+1:N+ps.M]
            # And the same for energies and ribosome fractions
            as = w*C[indT,(N+ps.M+1):(2*N+ps.M)] .+ (1-w)*C[indT-1,(N+ps.M+1):(2*N+ps.M)]
            ϕs = w*C[indT,(2*N+ps.M+1):end] .+ (1-w)*C[indT-1,(2*N+ps.M+1):end]
            # Then count the number of surviving species
            srv[i,j] = count(x->(x>0.0),pops)
            # Calculate rate of change in previous two steps of the dynamics
            dC = ((C[indT,1:N] .- C[indT-1,1:N])/(T[indT] - T[indT-1]))./(C[indT,1:N])
            # Two conditions, strain alive, strain not dying off
            a = pops .> 0.0
            b = dC .>= -5e-9
            # Count number of strains fufilling both conditions
            svinf[i,j] = count(x->x==2,a.+b)
            # Remove any negative metabolite concentrations
            for k = 1:ps.M
                if concs[k] < 0.0
                    concs[k] = 0.0
                end
            end
            # Now its safe to calculate the dissipation
            dsp[i,j] = dissipation(ps,ms,[pops;concs;as;ϕs])
            # Find dissipation rate per unit biomass
            b_dsp[i,j] = dsp[i,j]/sum(pops)
            # Now want to find and store all the fractional dissipations
            for k = 1:N
                # Need to determine the non-zero populations
                if pops[k] != 0.0
                    # Add population to vector of populations
                    pps[j] = cat(pps[j],log10(pops[k]),dims=1)
                    # Find vector of η values
                    ηs = ms[k].η
                    # make temporary vector for percentages
                    pdt = zeros(length(ηs))
                    # loop over η values
                    for l = 1:length(ηs)
                        # Find standard Gibbs free energy
                        dG = ps.reacs[ms[k].Reacs[l]].ΔG0
                        # Calculate percentage dissipated (under standard conditions)
                        pdt[l] = (ηs[l]*ΔGATP+dG)/(dG)
                    end
                    # Here I cat into the preallocated array of arrays
                    pd[j] = cat(pd[j],pdt,dims=1)
                    # Next want to catagorise this surviving microbes reactions
                    for l = 1:ms[k].R
                        # Find substrate concentration
                        S = concs[ps.reacs[ms[k].Reacs[l]].Rct]
                        # Check for case where the reaction has no substrate
                        if S == 0.0
                            inact[i,j] += 1.0
                        else
                            # Find product concentration
                            P = concs[ps.reacs[ms[k].Reacs[l]].Prd]
                            ΔG0 = ps.reacs[ms[k].Reacs[l]].ΔG0
                            # Find thermodynamic inhibition factor for this reaction
                            θ1 = θ(S,P,ps.T,ms[k].η[l],ΔG0)
                            # Use this factor to find obligate and non-obligate reactions
                            if θ1 >= 0.99
                                inact[i,j] += 1.0
                            else
                                # Find thermodynamic inhibition factor in case where product has built up
                                P = ps.κ[1]/ps.δ[1]
                                θ2 = θ(S,P,ps.T,ms[k].η[l],ΔG0)
                                if θ2 >= 0.5
                                    ob[i,j] += 1.0
                                else
                                    nob[i,j] += 1.0
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    # Set up plotting
    pyplot()
    # Set a color-blind friendly palette
    theme(:wong2,dpi=200)
    # Preallocate vector of labels
    Tlbs = Array{String,2}(undef,1,length(Ts))
    # Insert all elements to it
    for i = 1:length(Ts)
        Tlbs[i] = "T = $(Ts[i])s"
    end
    # Now plot the histograms
    m1 = L"^{-1}"
    # Pick appropriate indices
    inds = [3,5]
    histogram(srv[:,inds],labels=Tlbs[:,inds],fillalpha=0.75)
    plot!(title="$(R) reactions",xlabel="Number of surviving strains")
    savefig("Output/SvType$(R).png")
    inds = [5]
    histogram(svinf[:,inds],labels="",fillalpha=0.75)
    plot!(title="$(R) reactions",xlabel="Survivors at infinity")
    savefig("Output/InfSvType$(R).png")
    inds = [2,5]
    histogram(dsp[:,inds],labels=Tlbs[:,inds],fillalpha=0.75)
    plot!(title="$(R) reactions",xlabel="Dissipation rate (Js$(m1))")
    savefig("Output/DispType$(R).png")
    inds = [1,5]
    histogram(b_dsp[:,inds],labels=Tlbs[:,inds],fillalpha=0.75)
    plot!(title="$(R) reactions",xlabel="Dissipation rate (Js$(m1)) per cell")
    savefig("Output/BDispType$(R).png")
    inds = [1,5]
    histogram(pd[inds],labels=Tlbs[:,inds],fillalpha=0.75)
    plot!(title="$(R) reactions",xlabel="Percentage of free energy dissipated")
    savefig("Output/RDispType$(R).png")
    # Find and plot percentage active
    inds = [2,5]
    pa = 100*(1 .- inact[:,inds]./(inact[:,inds].+ob[:,inds].+nob[:,inds]))
    histogram(pa,labels=Tlbs[:,inds],fillalpha=0.75)
    plot!(title="$(R) reactions",xlabel="Percentage of reactions active")
    savefig("Output/PerActType$(R).png")
    po = 100*ob[:,inds]./(ob[:,inds].+nob[:,inds])
    histogram(po,labels=Tlbs[:,inds],fillalpha=0.75)
    plot!(title="$(R) reactions",xlabel="Percentage of reactions obligate")
    savefig("Output/PerObType$(R).png")
    # Final plot histogram of the populations (as measure of fitness)
    inds = [5]
    histogram(pps[inds],labels=Tlbs[:,inds],fillalpha=0.75)
    plot!(title="$(R) reactions",xlabel="Population (log scale)")
    savefig("Output/PopsType$(R).png")
    return(nothing)
end

# Function to analyse and plot the flux through networks
function net_vis()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 2
        error("need to specify community and number of repeats")
    end
    # Preallocate the variables I want to extract from the input
    R = 0
    nR = 0
    # Check that all arguments can be converted to integers
    try
        R = parse(Int64,ARGS[1])
        nR = parse(Int64,ARGS[2])
    catch e
           error("both inputs must be integer")
    end
    # Check that simulation type is valid
    if R < 1
        error("each strain must have more than 1 reaction")
    end
    # Check that number of simulations is greater than 0
    if nR < 1
        error("Number of repeats cannot be less than 1")
    end
    println("Compiled!")
    # Read in standard parameter file
    pfile = "Data/Type$(R)/ParasType$(R)Run1.jld"
    if ~isfile(pfile)
        error("run 1 is missing a parameter file")
    end
    ps = load(pfile,"ps")
    # Preallocate vectors to store data for histograms
    cA = zeros(nR)
    mf = zeros(nR)
    O = ps.O
    # Preallocate flows through reactions
    fR = zeros(ps.O)
    fRT = zeros(ps.O)
    # Loop over
    for i = 1:nR
        # Read in relevant files
        pfile = "Data/Type$(R)/ParasType$(R)Run$(i).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        ofile = "Data/Type$(R)/OutputType$(R)Run$(i).jld"
        if ~isfile(ofile)
            error("run $(i) is missing an output file")
        end
        efile = "Data/Type$(R)/ExtinctType$(R)Run$(i).jld"
        if ~isfile(efile)
            error("run $(i) is missing an extinct file")
        end
        # Basically just loading everything out as I'm not sure what I'll need
        ps = load(pfile,"ps")
        C = load(ofile,"C")
        T = load(ofile,"T")
        out = load(ofile,"out")
        ded = load(efile,"ded")
        # Loop over all reactions
        for j = 1:ps.O
            # Loop over microbes
            for k = 1:ps.N
                # Check if reaction exists in microbe
                if any(ps.mics[k].Reacs .== j) == true
                    # If so then find reaction
                    ind = findfirst(x->x==j,ps.mics[k].Reacs)
                    # Find final ribosome fraction for this strain
                    ϕR = out[2*ps.N+ps.M+k]
                    # Find amount of enzyme dedicated to reaction
                    E = Eα(ϕR,ps.mics[k],ind)
                    # Find substrate and product
                    S = out[ps.N+ps.reacs[j].Rct]
                    P = out[ps.N+ps.reacs[j].Prd]
                    # find reaction rate
                    q = qs(S,P,E,ind,ps.mics[k],ps.T,ps.reacs[j])
                    # Flux is reaction rate * population
                    fR[j] += q*out[k]
                end
            end
        end
        # Now want to put fluxes in units of moles
        fR = fR./NA
        fRT = fRT .+ fR
        # count the number of active reactions
        cA[i] = count(x->x>0.0,fR)
        # Find reaction with the maximum flux
        _, mf[i] = findmax(fR)
        # Plot bar charts for specfic individual ecosystems
        if i == 20
            # Set up plotting
            pyplot()
            # Set a color-blind friendly palette
            theme(:wong2,dpi=200)
            bar(fR,label="",xlabel="Reaction number",ylabel="Flux")
            plot!(title="$(R) reactions")
            savefig("Output/ReacsType$(R).png")
        end
    end
    # Set up plotting
    pyplot()
    # Set a color-blind friendly palette
    theme(:wong2,dpi=200)
    # Plot histogram of the number of active reactions
    histogram(cA,bins=range(0,stop=O+1,length=O+2))
    plot!(title="$(R) reactions",xlabel="Number of active reactions")
    savefig("Output/ActiveReactionsType$(R).png")
    histogram(mf,bins=range(1,stop=O+1,length=O+1))
    plot!(title="$(R) reactions",xlabel="Reaction with greatest flux")
    savefig("Output/MaxFluxType$(R).png")
    # Now find average flux
    fRT /= nR
    # Then plot
    bar(fRT,label="",xlabel="Reaction number",ylabel="Average flux")
    plot!(title="$(R) reactions")
    savefig("Output/AverageReacsType$(R).png")
    return(nothing)
end

# Function to make plots to show the shape of the thermodynamic tradeoff
function plt_trdff()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 2
        error("need to specify community and number of repeats")
    end
    # Preallocate the variables I want to extract from the input
    R = 0
    nR = 0
    # Check that all arguments can be converted to integers
    try
        R = parse(Int64,ARGS[1])
        nR = parse(Int64,ARGS[2])
    catch e
           error("both inputs must be integer")
    end
    # Check that simulation type is valid
    if R < 1
        error("each strain must have more than 1 reaction")
    end
    # Check that number of simulations is greater than 0
    if nR < 1
        error("Repeat number cannot be less than 1")
    end
    println("Compiled!")
    # Read in standard parameter file
    pfile = "Data/Type$(R)/ParasType$(R)Run$(nR).jld"
    if ~isfile(pfile)
        error("run $(nR) is missing a parameter file")
    end
    # Read in output data
    ofile = "Data/Type$(R)/OutputType$(R)Run$(nR).jld"
    if ~isfile(ofile)
        error("run $(nR) is missing an output file")
    end
    # Read in relevant data
    ps = load(pfile,"ps")
    out = load(ofile,"out")
    # Pick first reaction of the first microbe (guarenteed to exist)
    rc = ps.reacs[ps.mics[1].Reacs[1]]
    # Find substrate and product concentrations
    S = out[ps.N+rc.Rct]
    P = out[ps.N+rc.Prd]
    # Also store Gibbs free energy
    ΔG = rc.ΔG0
    # Amount of enzyme is chosen and fixed
    E = 100000.0
    # Set minimum equilibrium product to substrate ratio
    mratio = 1e-10
    # Use to generate appropriate set of η values
    ηmax = -(ΔG + Rgas*ps.T*log(mratio))/(ΔGATP)
    ηs = collect(range(0.33;length=200,stop=ηmax))
    # Preallocate vectors to store
    rs = zeros(length(ηs))
    θs = zeros(length(ηs))
    # Extra data to show impact of syntrophy
    Ps = [100.0*P,10.0*P,P,P/10.0,P/100.0]
    rs2 = zeros(length(ηs),length(Ps))
    # Loop over η values
    for i = 1:length(ηs)
        # Calculate thermodynamic inhibition
        θs[i] = θ_smooth(S,P,ps.T,ηs[i],ΔG)
        # Then use to calculate rate
        rs[i] = qs(ps.mics[1],S,P,E,θs[i])
        # Calculate the data needed to show the effect of syntrophy
        for j = 1:size(rs2,2)
            # Find temporary θ value for this case
            θt = θ_smooth(S,Ps[j],ps.T,ηs[i],ΔG)
            rs2[i,j] = qs(ps.mics[1],S,Ps[j],E,θt)
        end
    end
    # Make a vector of η*rate
    as = ηs.*rs
    # Rescale vectors as fractions of maximum
    rs = rs/maximum(rs)
    as2 = as/maximum(as)
    # Now setup plotting
    pyplot()
    # Set a color-blind friendly palette
    theme(:wong2,dpi=200)
    plot(ηs,rs,label="Reaction rate")
    plot!(ηs,θs,label="Inhibition")
    plot!(ηs,as2,label="ATP rate",xlabel=L"\eta")
    savefig("Output/TrdOff.png")
    # Do plot of just the tradeoff
    plot(ηs,as,label="",xlabel=L"\eta",ylabel="ATP production rate")
    savefig("Output/BareTrdOff.png")
    # Want a detailed visulisation of the peak
    inds = findall(x->(3.25<=x<=3.75),ηs)
    plot(ηs[inds],rs[inds],label="Reaction rate")
    plot!(ηs[inds],θs[inds],label="Inhibition")
    plot!(ηs[inds],as2[inds],label="ATP rate",xlabel=L"\eta")
    savefig("Output/PeakTrdOff.png")
    # Make a vector of η*rate
    as3 = ηs.*rs2
    # Define labels for the plot
    lbs = Array{String,2}(undef,1,5)
    for i = 1:5
        lbs[i] = "Prd conc = $(round(Ps[i],sigdigits=3))"
    end
    # Now calculate and plot syntrophy stuff
    plot(ηs,as3,xlabel=L"\eta",labels=lbs)
    savefig("Output/SynTrdOff.png")
    return(nothing)
end

# function to plot scatter graphs of effeciency verses expression
function react_scat()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 2
        error("need to specify community and number of repeats")
    end
    # Preallocate the variables I want to extract from the input
    R = 0
    nR = 0
    # Check that all arguments can be converted to integers
    try
        R = parse(Int64,ARGS[1])
        nR = parse(Int64,ARGS[2])
    catch e
           error("both inputs must be integer")
    end
    # Check that simulation type is valid
    if R < 1
        error("each strain must have more than 1 reaction")
    end
    # Check that number of simulations is greater than 0
    if nR < 1
        error("Number of repeats cannot be less than 1")
    end
    println("Compiled!")
    efT = Float64[]
    ϕpT = Float64[]
    actT = Bool[]
    # Now loop over repeats
    for i = 1:nR
        # First check that files exists
        pfile = "Data/Type$(R)/ParasType$(R)Run$(i).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        ofile = "Data/Type$(R)/OutputType$(R)Run$(i).jld"
        if ~isfile(ofile)
            error("run $(i) is missing an output file")
        end
        # Read in relevant data
        ps = load(pfile,"ps")
        out = load(ofile,"out")
        C = load(ofile,"C")
        T = load(ofile,"T")
        # Preallocate vector so that eventual extinct strains can be removed
        ext = BitArray(undef,ps.N)
        # Loop over every microbe
        for j = 1:ps.N
            # Find microbes index in full dynamics
            ind = findfirst(x->x==out[j],C[end,:])
            # Calculate rate change in final two steps of the dynamics
            dC = ((C[end,ind] - C[end-1,ind])/(T[end] - T[end-1]))/(C[end,ind])
            # If strain decreasing at greater than threshold rate then assume extinct
            if dC <= -5e-9
                ext[j] = 1
            else
                ext[j] = 0
            end
        end
        # Now need to find reduced parameter set
        ps = extinction(ps,ext)
        # Preallocate memory storage
        ef = zeros(ps.N*R) # Efficency
        ϕps = zeros(ps.N*R) # Relative activations
        act = fill(false,ps.N*R) # Reaction active
        # Loop over surviving strains
        for j = 1:ps.N
            # Then loop over each reaction
            for k = 1:R
                # Relative activation is easy to find
                ϕps[(j-1)*R+k] = ps.mics[j].ϕP[k]
                # First find standard Gibbs free energy
                dG = ps.reacs[ps.mics[j].Reacs[k]].ΔG0
                # Then use this to calculate the effeciency
                ef[(j-1)*R+k] = (ps.mics[j].η[k]*ΔGATP+dG)/(dG)
                # Find substrate and product indices
                indS = ps.reacs[ps.mics[j].Reacs[k]].Rct
                indP = ps.reacs[ps.mics[j].Reacs[k]].Prd
                # Find substrate and product concentration
                S = out[ps.N+sum(ext)+indS]
                P = out[ps.N+sum(ext)+indP]
                # Find thermodynamic inhibition factor for this reaction
                θ1 = θ(S,P,ps.T,ps.mics[j].η[k],dG)
                # Check the reaction has substrate and is thermodynamically feasible
                if S != 0.0 && θ1 <= 0.99
                    act[(j-1)*R+k] = true
                end
            end
        end
        # No loop is done cat results onto overall vector
        efT = cat(efT,ef,dims=1)
        ϕpT = cat(ϕpT,ϕps,dims=1)
        actT = cat(actT,act,dims=1)
    end
    # Now setup plotting
    pyplot()
    # Set a color-blind friendly palette
    theme(:wong2,dpi=200)
    scatter(efT[actT],ϕpT[actT],label="Active",xlabel="Efficency",ylabel="Activation",color=2)
    scatter!(efT[.~actT],ϕpT[.~actT],label="Inactive",color=3)
    savefig("Output/ScatterType$(R).png")
    # Count the number of reactions in each catgory
    c0 = count(x->(x<0.0),efT)
    c1 = count(x->(0.2>x>=0.0),efT)
    c2 = count(x->(0.4>x>=0.2),efT)
    c3 = count(x->(0.6>x>=0.4),efT)
    c4 = count(x->(0.8>x>=0.6),efT)
    c5 = count(x->(x>=0.8),efT)
    println("Super efficent = $(100*c0/length(efT))%")
    println("Highly efficent = $(100*c1/length(efT))%")
    println("Moderately efficent = $(100*c2/length(efT))%")
    println("Moderately inefficent = $(100*c3/length(efT))%")
    println("Highly inefficent = $(100*c4/length(efT))%")
    println("Super inefficent = $(100*c5/length(efT))%")
    return(nothing)
end

@time hist_time()
