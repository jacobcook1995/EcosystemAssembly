# This script exists to plot the output of form_comm.jl
using Assembly
using StatsPlots
using LaTeXStrings
using JLD
using StatsBase
using Statistics
using DataFrames
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

# A hardcoded extended uniform distribution function used by net_vis()
function f_ext_un(x::StepRangeLen,syn::Bool)
    # Preallocate output
    y = zeros(length(x))
    # Change limits based on whether syntrophy is on or not
    if syn == true
        h1 = 1.0
        h2 = 1.0
    else
        h1 = 0.7352
        h2 = 0.8624
    end
    for i = 1:length(x)
        # Check if it's high enough to be in first uniform range
        if x[i] <= h1
            y[i] += 1
        end
        # Check if it's high enough to be in second uniform range
        if x[i] <= h2
            y[i] += 1
        end
        # Then subtract if too low
        if x[i] < 0.08756
            y[i] -= 1
        end
        if x[i] < 0.1685
            y[i] -= 1
        end
    end
    return(y)
end

# Function to analyse and plot the flux through networks
function net_vis()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 4
        error("Insufficent inputs provided (looking for 4)")
    end
    # Preallocate the variables I want to extract from the input
    Rl = 0
    Ru = 0
    syn = true
    nR = 0
    # Check that all arguments can be converted to integers
    try
        Rl = parse(Int64,ARGS[1])
        Ru = parse(Int64,ARGS[2])
        syn = parse(Bool,ARGS[3])
        nR = parse(Int64,ARGS[4])
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
    if nR < 1
        error("simulation number can't be less than 1")
    end
    println("Compiled!")
    # Read in standard parameter file
    pfile = "Data/$(Rl)-$(Ru)$(syn)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run1.jld"
    if ~isfile(pfile)
        error("run 1 is missing a parameter file")
    end
    ps = load(pfile,"ps")
    # Preallocate vectors to store data for histograms
    cA = zeros(Int64,nR)
    mf = zeros(Int64,nR)
    O = ps.O
    # Preallocate flows through reactions
    fRT = zeros(ps.O)
    # Mass renormalised version
    fRTm = zeros(ps.O)
    # Percentage flow through most contributing strain
    prf = zeros(nR,ps.O)
    # Preallocate storage for best kinetic parameters
    KSs = []
    krs = []
    kcs = []
    ϕps = []
    efs = []
    # And for the full kinetic parameters
    KSF = []
    krF = []
    kcF = []
    ϕpF = []
    efF = []
    # mimimum product to substrate ratio (to calculate) the efficency
    mr = 1e-2
    # To store number of ecosystems that reaction is present in
    pres = zeros(ps.O)
    # Data structure to store abundances
    abnd = zeros(ps.O,1)
    # Loop over repeats
    for i = 1:nR
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
        # Basically just loading everything out as I'm not sure what I'll need
        ps = load(pfile,"ps")
        C = load(ofile,"C")
        T = load(ofile,"T")
        out = load(ofile,"out")
        ded = load(efile,"ded")
        # Temporary way of storing the fluxes for this realisation
        fR = zeros(ps.O)
        # Mass renormalised version
        fRm = zeros(ps.O)
        # Highest flux from one reaction (for particular strain)
        bstf = zeros(ps.N)
        # reaction that acheives this
        bstr = zeros(Int64,ps.N)
        # Highest flux through one strain
        hghf = zeros(ps.O)
        # microbes that acheive that flux
        bstm = zeros(Int64,ps.O)
        # Quickly loop over microbes to find and store kinetics data
        for j = 1:ps.N
            KSF = cat(KSF,ps.mics[j].KS,dims=1)
            kcF = cat(kcF,ps.mics[j].kc,dims=1)
            krF = cat(krF,ps.mics[j].kr,dims=1)
            ϕpF = cat(ϕpF,ps.mics[j].ϕP,dims=1)
            # Find vector of ΔG0 values
            dG = ps.reacs[ps.mics[j].Reacs].↦:ΔG0
            # Use to calculate percentage retained (under standard conditions)
            efT = -(ps.mics[j].η*ΔGATP)./(dG.+Rgas*ps.T*log(mr))
            # cat efficency in
            efF = cat(efF,efT,dims=1)
        end
        # Make temp vectors for kinetic parameters
        ϕp = zeros(ps.N)
        KS = zeros(ps.N)
        kr = zeros(ps.N)
        kc = zeros(ps.N)
        ef = zeros(ps.N)
        # Loop over all reactions
        for j = 1:ps.O
            # Check if reaction is present
            if any(j .∈ (ps.mics.↦:Reacs))
                pres[j] += 1
                # Count the number of occurances
                n = count(j .∈ (ps.mics.↦:Reacs))
                # Update container until its the right size
                while n + 1 > size(abnd,2)
                    tmp = zeros(ps.O)
                    abnd = cat(abnd,tmp,dims=2)
                end
                # Add to correct abundance
                abnd[j,n+1] += 1
            else
                abnd[j,1] += 1
            end
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
                    # Mass renormalised version
                    fRm[j] += q
                    # Find eta value for strain
                    ηf = ps.mics[k].η[ind]
                    # Check if ATP flux is higher than other reactions for strain
                    if ηf*q*out[k] > bstf[k]
                        bstf[k] = ηf*q*out[k]
                        bstr[k] = j
                    end
                    # Check if reaction flux is higher than for other strains
                    if q*out[k] > hghf[j]
                        hghf[j] = q*out[k]
                        bstm[j] = k
                    end
                end
            end
        end
        # Now want to find percentage of total flux carried by greatest contributing strain
        for j = 1:ps.O
            if fR[j] > 0.0
                prf[i,j] = (hghf[j]/fR[j])*100.0
            end
        end
        # Now for each strain want to find relevant parameters
        for j = 1:ps.N
            # Find where in the vector the reaction number is found
            ind = findfirst(x->(x==bstr[j]),ps.mics[j].Reacs)
            KS[j] = ps.mics[j].KS[ind]
            kc[j] = ps.mics[j].kc[ind]
            kr[j] = ps.mics[j].kr[ind]
            ϕp[j] = ps.mics[j].ϕP[ind]
            # Find relevant ΔG0 value
            dG = ps.reacs[ps.mics[j].Reacs[ind]].ΔG0
            # Use to calculate percentage dissipated (under standard conditions)
            ef[j] = -(ps.mics[j].η[ind]*ΔGATP)/(dG+Rgas*ps.T*log(mr))
        end
        # Add these kinetic parameters to storage
        KSs = cat(KSs,KS,dims=1)
        kcs = cat(kcs,kc,dims=1)
        krs = cat(krs,kr,dims=1)
        ϕps = cat(ϕps,ϕp,dims=1)
        efs = cat(efs,ef,dims=1)
        # Now want to put fluxes in units of moles
        fR = fR./NA
        fRm = fRm
        # Then add to the total
        fRT = fRT .+ fR
        fRTm = fRTm .+ fRm
        # count the number of active reactions
        cA[i] = count(x->x>0.0,fR)
        # Find reaction with the maximum flux
        _, mf[i] = findmax(fR)
    end
    # Set up plotting
    pyplot()
    # Set a color-blind friendly palette
    theme(:wong2,dpi=200)
    # Make tick labels
    xs = Array{String,1}(undef,length(fRT))
    for i = 1:length(xs)
        # Find substrate and product numbers
        s = ceil(Int64,i/2)
        p = 2 + floor(Int64,i/2)
        xs[i] = "$s→$(p)"
    end
    tl = ""
    if syn == true
        tl = "$(Rl)-$(Ru) reactions per strain"
    else
        tl = "$(Rl)-$(Ru) reactions per strain (no syntrophy)"
    end
    # Plot "presence" data
    groupedbar([nR.-pres pres],bar_position=:stack,labels=["Absent" "Present"])
    plot!(xticks=(1:ps.O,xs),xlabel="Reaction",ylabel="Number of ecosystems")
    plot!(title=tl,legend=:outerright)
    savefig("Output/$(Rl)-$(Ru)$(syn)/PresenseType$(Rl)-$(Ru)$(syn).png")
    # Now plot "abundance" data
    lbs = Array{String,2}(undef,1,size(abnd,2))
    for i = 1:length(lbs)
        lbs[i] = string(i-1)
    end
    groupedbar(abnd,bar_position=:stack,labels=lbs,legend=:outerright)
    plot!(xticks=(1:ps.O,xs),xlabel="Reaction",ylabel="Number of ecosystems")
    plot!(title=tl)
    savefig("Output/$(Rl)-$(Ru)$(syn)/AbundanceType$(Rl)-$(Ru)$(syn).png")
    # Plot histogram of the number of active reactions
    histogram(cA,bins=range(0,stop=O+1,length=O+2))
    plot!(title=tl,xlabel="Number of active reactions")
    savefig("Output/$(Rl)-$(Ru)$(syn)/ActiveReactionsType$(Rl)-$(Ru)$(syn).png")
    histogram(mf,bins=range(1,stop=O+1,length=O+1),xticks=(1:ps.O,xs))
    plot!(title=tl,xlabel="Reaction with greatest flux")
    savefig("Output/$(Rl)-$(Ru)$(syn)/MaxFluxType$(Rl)-$(Ru)$(syn).png")
    # Now find average flux
    fRT /= nR
    # Then plot
    bar(fRT,xticks=(1:ps.O,xs),label="",xlabel="Reaction",ylabel="Average flux")
    plot!(title=tl)
    savefig("Output/$(Rl)-$(Ru)$(syn)/AverageReacsType$(Rl)-$(Ru)$(syn).png")
    # Then find and plot average mass renormalised flux
    fRTm /= nR
    bar(fRTm,xticks=(1:ps.O,xs),label="",xlabel="Reaction",ylabel="Mass renormalised average flux")
    plot!(title=tl)
    savefig("Output/$(Rl)-$(Ru)$(syn)/AverageReacsMassType$(Rl)-$(Ru)$(syn).png")
    # Preallocate percent data to plot
    perc = zeros(ps.O)
    perc2 = zeros(ps.O)
    for i = 1:ps.O
        # Count zero elements and remove from sum
        ct = count(x->x==0.0,prf[:,i])
        if ct != nR
            perc[i] = sum(prf[:,i])/(nR-ct)
        else
            perc[i] = 0.0
        end
        # Count single elements and remove from sum
        ct2 = count(x->x==100.0,prf[:,i])
        if ct2 != 0
            perc2[i] = (sum(prf[:,i])-100.0*ct2)/(nR-ct-ct2)
        else
            perc2[i] = sum(prf[:,i])/(nR-ct)
        end
    end
    # Now plot how many strains have each reaction
    bar(perc,label="",ylabel="Average percentage of flux through top strain")
    plot!(title=tl,xlabel="Reaction",xticks=(1:ps.O,xs))
    savefig("Output/$(Rl)-$(Ru)$(syn)/PercentFluxType$(Rl)-$(Ru)$(syn).png")
    bar(perc2,label="",ylabel="Percentange of flux (ignoring singles)")
    plot!(title=tl,xlabel="Reaction",xticks=(1:ps.O,xs))
    savefig("Output/$(Rl)-$(Ru)$(syn)/RedPercentFluxType$(Rl)-$(Ru)$(syn).png")
    # Define distribution functions
    @. f_ln(x,μ,σ) = (1/(x*σ*sqrt(2*π)))*exp(-((log(x)-μ)^2)/(2*σ^2))
    @. f_nm(x,μ,σ) = (1/(σ*sqrt(2*π)))*exp(-(1/2)*((x-μ)/(σ))^2)
    # Then manually calculate histograms
    bins1 = range(1e-4,stop=maximum(KSs),length=500)
    h1 = fit(Histogram,KSs,bins1,closed=:right)
    bins2 = range(1e-4,stop=maximum(kcs),length=500)
    h2 = fit(Histogram,kcs,bins2,closed=:right)
    bins3 = range(1e-4,stop=maximum(krs),length=500)
    h3 = fit(Histogram,krs,bins3,closed=:right)
    bins4 = range(1e-4,stop=maximum(ϕps),length=500)
    h4 = fit(Histogram,ϕps,bins4,closed=:right)
    bins5 = range(-0.0,stop=1.00,length=500)
    h5 = fit(Histogram,efs,bins5,closed=:right)
    # Then plot as bar charts, with inital distribution included
    bar(h1,label="Main reaction")
    # Find variables needed for the distribution
    μ1 = log((1/4)*5.5e-3)
    σ1 = log(2)
    # Renormalise distribution for plotting
    d1 = f_ln(bins1,μ1,σ1)
    d1 /= maximum(d1)
    plot!(bins1,maximum(h1.weights)*d1,label="Orginal distribution")
    plot!(title=tl,xlabel=L"K_S")
    savefig("Output/$(Rl)-$(Ru)$(syn)/TopKSType$(Rl)-$(Ru)$(syn).png")
    bar(h2,label="Main reaction")
    μ2 = log(10.0)
    σ2 = log(2)
    d2 = f_ln(bins2,μ2,σ2)
    d2 /= maximum(d2)
    plot!(bins2,maximum(h2.weights)*d2,label="Orginal distribution")
    plot!(title=tl,xlabel=L"k_c")
    savefig("Output/$(Rl)-$(Ru)$(syn)/TopkcType$(Rl)-$(Ru)$(syn).png")
    bar(h3,label="Main reaction")
    μ3 = log(10.0)
    σ3 = log(2)
    d3 = f_ln(bins3,μ3,σ3)
    d3 /= maximum(d3)
    plot!(bins3,maximum(h3.weights)*d3,label="Orginal distribution")
    plot!(title=tl,xlabel=L"k_r")
    savefig("Output/$(Rl)-$(Ru)$(syn)/TopkrType$(Rl)-$(Ru)$(syn).png")
    bar(h4,label="Main reaction")
    plot!(title=tl,xlabel=L"\phi_p")
    savefig("Output/$(Rl)-$(Ru)$(syn)/TopPhiPType$(Rl)-$(Ru)$(syn).png")
    bar(h5,label="Main reaction")
    d5 = f_ext_un(bins5,syn)
    d5 /= maximum(d5)
    plot!(bins5,maximum(h5.weights)*d5,label="Orginal distribution",legend=:topleft)
    plot!(title=tl,xlabel="Fraction of free energy retained")
    savefig("Output/$(Rl)-$(Ru)$(syn)/TopEfficencyType$(Rl)-$(Ru)$(syn).png")
    # Repeat process for full kinetics data
    bins1 = range(1e-4,stop=maximum(KSF),length=500)
    h1 = fit(Histogram,KSF,bins1,closed=:right)
    bins2 = range(1e-4,stop=maximum(kcF),length=500)
    h2 = fit(Histogram,kcF,bins2,closed=:right)
    bins3 = range(1e-4,stop=maximum(krF),length=500)
    h3 = fit(Histogram,krF,bins3,closed=:right)
    bins4 = range(1e-4,stop=maximum(ϕpF),length=500)
    h4 = fit(Histogram,ϕpF,bins4,closed=:right)
    bins5 = range(0.0,stop=1.00,length=500)
    h5 = fit(Histogram,efF,bins5,closed=:right)
    # Then plot as bar charts, with inital distribution included
    bar(h1,label="Main reaction")
    # Find variables needed for the distribution
    μ1 = log((1/4)*5.5e-3)
    σ1 = log(2)
    # Renormalise distribution for plotting
    d1 = f_ln(bins1,μ1,σ1)
    d1 /= maximum(d1)
    plot!(bins1,maximum(h1.weights)*d1,label="Orginal distribution")
    plot!(title=tl,xlabel=L"K_S")
    savefig("Output/$(Rl)-$(Ru)$(syn)/AllKSType$(Rl)-$(Ru)$(syn).png")
    bar(h2,label="All reactions")
    μ2 = log(10.0)
    σ2 = log(2)
    d2 = f_ln(bins2,μ2,σ2)
    d2 /= maximum(d2)
    plot!(bins2,maximum(h2.weights)*d2,label="Orginal distribution")
    plot!(title=tl,xlabel=L"k_c")
    savefig("Output/$(Rl)-$(Ru)$(syn)/AllkcType$(Rl)-$(Ru)$(syn).png")
    bar(h3,label="All reactions")
    μ3 = log(10.0)
    σ3 = log(2)
    d3 = f_ln(bins3,μ3,σ3)
    d3 /= maximum(d3)
    plot!(bins3,maximum(h3.weights)*d3,label="Orginal distribution")
    plot!(title=tl,xlabel=L"k_r")
    savefig("Output/$(Rl)-$(Ru)$(syn)/AllkrType$(Rl)-$(Ru)$(syn).png")
    bar(h4,label="All reactions")
    plot!(title=tl,xlabel=L"\phi_p")
    savefig("Output/$(Rl)-$(Ru)$(syn)/AllPhiPType$(Rl)-$(Ru)$(syn).png")
    bar(h5,label="All reactions")
    d5 = f_ext_un(bins5,syn)
    d5 /= maximum(d5)
    plot!(bins5,maximum(h5.weights)*d5,label="Orginal distribution",legend=:topleft)
    plot!(title=tl,xlabel="Fraction of free energy retained")
    savefig("Output/$(Rl)-$(Ru)$(syn)/AllEfficencyType$(Rl)-$(Ru)$(syn).png")
    return(nothing)
end

# Function to make plots to show the shape of the thermodynamic tradeoff
function plt_trdff()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 4
        error("Insufficent inputs provided (looking for 4)")
    end
    # Preallocate the variables I want to extract from the input
    Rl = 0
    Ru = 0
    syn = true
    nR = 0
    # Check that all arguments can be converted to integers
    try
        Rl = parse(Int64,ARGS[1])
        Ru = parse(Int64,ARGS[2])
        syn = parse(Bool,ARGS[3])
        nR = parse(Int64,ARGS[4])
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
    if nR < 1
        error("simulation number can't be less than 1")
    end
    println("Compiled!")
    # Read in relevant files
    pfile = "Data/$(Rl)-$(Ru)$(syn)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(nR).jld"
    if ~isfile(pfile)
        error("run $(nR) is missing a parameter file")
    end
    ofile = "Data/$(Rl)-$(Ru)$(syn)/RedOutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(nR).jld"
    if ~isfile(ofile)
        error("run $(nR) is missing an output file")
    end
    efile = "Data/$(Rl)-$(Ru)$(syn)/RedExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(nR).jld"
    if ~isfile(efile)
        error("run $(nR) is missing an extinct file")
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
    wongc = get_color_palette(wong_palette,57)
    plot(ηs,rs,label="Reaction rate")
    plot!(ηs,θs,label="Inhibition")
    plot!(ηs,as2,label="ATP rate",xlabel="ATP per reaction event")
    savefig("Output/TrdOff.png")
    # Do plot of just the tradeoff
    plot(ηs,as,label="",xlabel="ATP per reaction event",ylabel="ATP production rate")
    savefig("Output/BareTrdOff.png")
    # Want a detailed visulisation of the peak
    inds = findall(x->(3.25<=x<=3.75),ηs)
    plot(ηs[inds],rs[inds],label="Reaction rate")
    plot!(ηs[inds],θs[inds],label="Inhibition")
    plot!(ηs[inds],as2[inds],label="ATP rate",xlabel="ATP per reaction event")
    savefig("Output/PeakTrdOff.png")
    # Make a vector of η*rate
    as3 = ηs.*rs2
    # Define labels for the plot
    lbs = Array{String,2}(undef,1,5)
    for i = 1:5
        lbs[i] = "Prd conc = $(round(Ps[i],sigdigits=3))"
    end
    # Now calculate and plot syntrophy stuff
    plot(ηs,as3,xlabel="ATP per reaction event",ylabel="ATP production rate",labels="")#lbs)
    savefig("Output/SynTrdOff.png")
    plot(ηs,θs,yscale=:log10,xlabel=L"\eta",ylabel=L"\theta",label="",ylims=(1e-15,1e1),color=wongc[2])
    savefig("Output/LogInhib.png")
    return(nothing)
end

# function I can use to test particular parameter sets
function test()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 6
        error("Insufficent inputs provided (looking for 6)")
    end
    # Preallocate the variables I want to extract from the input
    Rl = 0
    Ru = 0
    syn = true
    Nr = 0
    Ns = 0
    en = ARGS[6]
    # Check that all arguments can be converted to integers
    try
        Rl = parse(Int64,ARGS[1])
        Ru = parse(Int64,ARGS[2])
        syn = parse(Bool,ARGS[3])
        Nr = parse(Int64,ARGS[4])
        Ns = parse(Int64,ARGS[5])
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
    if Nr < 1
        error("simulation number can't be less than 1")
    end
    if Ns < 1
        error("number of strains can't be less than 1")
    end
    println("Compiled!")
    # Read in relevant files
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
    ps = load(pfile,"ps")
    C = load(ofile,"C")
    T = load(ofile,"T")
    out = load(ofile,"out")
    ded = load(efile,"ded")
    # ind = findfirst(x->x==out[1],C[end,:])
    # m = ps.mics[1]
    # Now move onto plotting
    pyplot()
    theme(:wong2,dpi=200)
    plot(yaxis=:log10)
    # Find indicies
    is = zeros(Int64,ps.N)
    for i = 1:ps.N
        is[i] = findfirst(x->x==out[i],C[end,:])
    end
    # Find orginal number of strains
    N = ps.N + length(ded)
    # Plot populations of final survivors
    plot(title="Final populations",yaxis=:log10)
    for i = is
         # Find and eliminate zeros so that they can be plotted on a log plot
         inds = (C[:,i] .> 0)
         plot!(T[inds],C[inds,i],label="")
    end
    savefig("Output/redpops.png")
    plot(T,C[:,(N+ps.M.+is)],label="")
    savefig("Output/redATPs.png")
    plot(T,C[:,(2*N+ps.M.+is)],label="")
    savefig("Output/redfracs.png")
    # Plot all the populations
    plot(title="Full populations",yaxis=:log10)
    for i = 1:N
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0)
        plot!(T[inds],C[inds,i],label="")
    end
    savefig("Output/fullpops.png")
    plot(T,C[:,(N+1):(N+ps.M)],label="")
    savefig("Output/concs.png")
    plot(T,C[:,(N+ps.M+1):(2*N+ps.M)],label="")
    savefig("Output/fullATPs.png")
    plot(T,C[:,(2*N+ps.M+1):(3*N+ps.M)],label="")
    savefig("Output/fullfracs.png")
    return(nothing)
end

# Function to plot distribution of fluxes and abundances
function flux_abund()
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
    # Data structure to store abundances
    abnd = []
    # Data structures to store fluxes
    fl = []
    flm = []
    Afl = []
    Aflm = []
    # Loop over repeats
    for i = 1:nR
        # Read in relevant files
        pfile = "Data/Type$(R)/RedParasType$(R)Run$(i).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        ofile = "Data/Type$(R)/RedOutputType$(R)Run$(i).jld"
        if ~isfile(ofile)
            error("run $(i) is missing an output file")
        end
        efile = "Data/Type$(R)/RedExtinctType$(R)Run$(i).jld"
        if ~isfile(efile)
            error("run $(i) is missing an extinct file")
        end
        # Basically just loading everything out as I'm not sure what I'll need
        ps = load(pfile,"ps")
        C = load(ofile,"C")
        T = load(ofile,"T")
        out = load(ofile,"out")
        ded = load(efile,"ded")
        # Add abundances to the vector
        abnd = cat(abnd,out[1:ps.N],dims=1)
        # Preallocate fluxes
        f = zeros(R*ps.N)
        fm = zeros(R*ps.N)
        Af = zeros(R*ps.N)
        Afm = zeros(R*ps.N)
        # Loop over all reactions
        for k = 1:ps.N
            # Loop over microbes
            c = 0
            for j = 1:ps.O
                # Check if reaction exists in microbe
                if any(ps.mics[k].Reacs .== j) == true
                    # Increment counter
                    c += 1
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
                    f[(k-1)*R+c] = q*out[k]
                    # Mass renormalised version
                    fm[(k-1)*R+c] = q
                    # Find eta value for strain
                    ηf = ps.mics[k].η[ind]
                    # ATP fluxes
                    Af[(k-1)*R+c] = q*out[k]*ηf
                    Afm[(k-1)*R+c] = q*ηf
                end
            end
        end
        # Cat vector of fluxes in
        fl = cat(fl,f,dims=1)
        flm = cat(flm,fm,dims=1)
        Afl = cat(Afl,Af,dims=1)
        Aflm = cat(Aflm,Afm,dims=1)
    end
    # Set up plotting
    pyplot()
    # Set a color-blind friendly palette
    theme(:wong2,dpi=200)
    histogram(abnd,label="",title="$(R) reactions per strain",xlabel="Abundances")
    savefig("Output/Type$(R)/StrainAbundType$(R).png")
    histogram(log10.(abnd),label="",title="$(R) reactions per strain",xlabel="Log Abundances")
    savefig("Output/Type$(R)/LogStrainAbundType$(R).png")
    histogram(fl,label="",title="$(R) reactions per strain",xlabel="Flux")
    savefig("Output/Type$(R)/FluxesType$(R).png")
    histogram(log10.(fl),label="",title="$(R) reactions per strain",xlabel="Log Flux")
    savefig("Output/Type$(R)/LogFluxesType$(R).png")
    histogram(flm,label="",title="$(R) reactions per strain",xlabel="Mass specific flux")
    savefig("Output/Type$(R)/FluxesMType$(R).png")
    histogram(log10.(flm),label="",title="$(R) reactions per strain",xlabel="Log Mass specific flux")
    savefig("Output/Type$(R)/LogFluxesMType$(R).png")
    histogram(Afl,label="",title="$(R) reactions per strain",xlabel="ATP flux")
    savefig("Output/Type$(R)/ATPFType$(R).png")
    histogram(log10.(Afl),label="",title="$(R) reactions per strain",xlabel="Log ATP flux")
    savefig("Output/Type$(R)/LogATPFType$(R).png")
    histogram(Aflm,label="",title="$(R) reactions per strain",xlabel="Mass specific ATP flux")
    savefig("Output/Type$(R)/ATPFMType$(R).png")
    histogram(log10.(Aflm),label="",title="$(R) reactions per strain",xlabel="Log Mass specific ATP flux")
    savefig("Output/Type$(R)/LogATPFMType$(R).png")
    return(nothing)
end

# Function to plot the basic information about the simulations
function basic_info()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 6
        error("Insufficent inputs provided (looking for 6)")
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
        Rl = parse(Int64,ARGS[1])
        Ru = parse(Int64,ARGS[2])
        syn = parse(Bool,ARGS[3])
        rps = parse(Int64,ARGS[4])
        Ni = parse(Int64,ARGS[5])
    catch e
        error("need to provide 4 integers, and a bool")
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
        error("number of initial strains can't be less than 1")
    end
    println("Compiled!")
    # Preallocate data
    svs = zeros(Int64,rps) # Number of survivors
    mbs = zeros(Int64,rps) # Number of metabolites
    hmb = zeros(Int64,rps) # Lowest (in energy hierachy) metabolite
    # Make containers for the data of uncertain size
    rcs = Array{Int64,1}(undef,0) # Number of reactions (per strain)
    Sbs = Array{Array{Int64,1},1}(undef,Ru-Rl+1) # Substrates used in reactions
    abds = [] # Abundances
    rabs = [] # Relative abundances
    effs = [] # Reaction efficencies
    # Defining stuff that then gets defined in the loop
    fnd = []
    avR = []
    # Number of metabolites
    M = 0
    # mimimum product to substrate ratio (to calculate) the efficency
    mr = 1e-2
    # Bools to check if data has been written in for each reaction number
    rs2 = fill(false,Ru-Rl+1)
    # Loop over repeats
    for i = 1:rps
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
        ps = load(pfile,"ps")
        out = load(ofile,"out")
        inf_out = load(ofile,"inf_out")
        # Calculate mean reaction value for ecosystem
        mR = sum(ps.mics.↦:R)/ps.N
        # Save metabolite number from the first parameter set
        if i == 1
            M = ps.M
            # Use to preallocate structure to store average R's
            avR = Array{Array{Float64,1},1}(undef,ps.M-1)
            # Vector to record if metabolite has previously been found
            fnd = fill(false,ps.M-1)
        end
        # To store which metabolites exist
        exst = fill(false,ps.M-1)
        # Check if each metabolite exists
        for j = 1:ps.M-1
            if inf_out[ps.N+j] > 0.0
                exst[j] = true
            end
        end
        # Count number of existing substrates
        ind = count(x->x==true,exst)
        # Store mean value of this index
        if fnd[ind] == false
            avR[ind] = [mR]
            fnd[ind] = true
        else
            avR[ind] = cat(avR[ind],mR,dims=1)
        end
        # Find and store number of survivors
        svs[i] = ps.N
        # Find and store number of reactions
        if ps.N > 0
            rcs = cat(rcs,ps.mics.↦:R,dims=1)
            # Check if data has already been saved
            rs = fill(false,Ru-Rl+1)
            # Make temporary vector to store substrates
            sbt = Array{Array{Int64,1},1}(undef,Ru-Rl+1)
            # Fill out temporary vector
            for j = 1:ps.N
                # Find reaction number for strain
                R = ps.mics[j].R - Rl + 1
                # Add to relevant collection, save for first time
                if rs[R] == false
                    sbt[R] = ps.reacs[ps.mics[j].Reacs].↦:Rct
                    rs[R] = true
                else
                    # Then cat for later times
                    sbt[R] = cat(sbt[R],ps.reacs[ps.mics[j].Reacs].↦:Rct,dims=1)
                end
            end
            # And cat into larger collection
            for j = 1:(Ru-Rl+1)
                # Similar step to ensure that first step is written in
                if rs2[j] == false && rs[j] == true
                    Sbs[j] = sbt[j]
                    rs2[j] = true
                elseif rs[j] == true
                    Sbs[j] = cat(Sbs[j],sbt[j],dims=1)
                end
            end
        end
        # Find and store abundances
        abds = cat(abds,inf_out[1:ps.N],dims=1)
        # And relative abundances
        rabs = cat(rabs,log10.(inf_out[1:ps.N]./sum(inf_out[1:ps.N])),dims=1)
        # Loop over microbes, to find efficencies
        for j = 1:ps.N
            # Find vector of ΔG0 values
            dG = ps.reacs[ps.mics[j].Reacs].↦:ΔG0
            # Use to calculate percentage dissipated (under standard conditions)
            efT = -(ps.mics[j].η.*ΔGATP)./(dG.+Rgas*ps.T*log(mr))
            # cat efficency in
            effs = cat(effs,efT,dims=1)
        end
        # Loop over metabolites to find those with non-zero concentrations
        cm = 0 # Set up counter
        mm = 0 # Lowest metabolite
        for j = 1:ps.M
            if inf_out[ps.N+j] > 0.0
                # Increment counter
                cm += 1
                # And save new minimum metabolite
                mm = j
            end
        end
        # Save results to vector
        mbs[i] = cm
        hmb[i] = mm
    end
    println("Extracted Data!")
    # Preallocate vectors
    ms = zeros(Ru-Rl+1)
    sds = zeros(Ru-Rl+1)
    Rs = collect(Rl:Ru)
    # Loop over number of reactions
    for j = Rl:Ru
        # Find indices of this reactions
        inds = findall(x->x==j,rcs)
        # Use to find mean and standard deviation of relevant abundances
        ms[j+1-Rl] = mean(abds[inds])
        sds[j+1-Rl] = std(abds[inds])
    end
    # Set up plotting
    pyplot()
    theme(:wong2,dpi=200)
    # Preallocate title
    tl = ""
    # Then choose title
    if syn == true
        tl = "$(Rl)-$(Ru) reactions per strain"
    else
        tl = "$(Rl)-$(Ru) reactions per strain (no syntrophy)"
    end
    # Make survivor bins
    bs = range(-0.25,stop=M-0.75,length=2*M)
    # And efficency bins
    ebs = range(0.0,stop=100.0,length=100)
    # And metabolite bins
    mbn = range(-0.25,stop=M+0.25,length=2*(M+1))
    # Add reaction bins
    rbs = range(Rl-0.25,stop=Ru+0.25,length=2*(Ru-Rl+1))
    # Plot histograms of the data
    histogram(svs,bins=bs,label="",xlabel="Number of strains",title=tl)
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/Survivors$(Rl)-$(Ru)$(syn)$(Ni).png")
    histogram(rcs,bins=rbs,label="",xlabel="Number of reactions",title=tl)
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/Reactions$(Rl)-$(Ru)$(syn)$(Ni).png")
    histogram(log10.(abds),label="",xlabel="Species abundance (log of number of cells)",title=tl)
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/Abundance$(Rl)-$(Ru)$(syn)$(Ni).png")
    # Rescale distribution to mean zero
    rrabs = rabs .- mean(rabs)
    # Then plot rescaled abundances
    histogram(rrabs,label="",xlabel="Rescaled log relative abundances",title=tl)
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RelAbundance$(Rl)-$(Ru)$(syn)$(Ni).png")
    histogram(effs*100.0,bins=ebs,label="",xlabel="Efficency",title=tl)
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/Efficency$(Rl)-$(Ru)$(syn)$(Ni).png")
    scatter([Rs],[ms],yerror=sds,label="")
    plot!(title=tl,xlabel="Number of reactions",ylabel="Strain abundance (number of cells)")
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/AbvsRct$(Rl)-$(Ru)$(syn)$(Ni).png")
    histogram(mbs,bins=mbn,label="",xlabel="Final number of metabolites",title=tl)
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/NoMets$(Rl)-$(Ru)$(syn)$(Ni).png")
    histogram(hmb,bins=mbn,label="",xlabel="Lowest energy metabolite",title=tl)
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/LowMet$(Rl)-$(Ru)$(syn)$(Ni).png")
    # Make complete collection of substrates
    SbsT = Array{Int64,1}(undef,0)
    # Plot number of reactions using each substrate for each reaction number
    for j = 1:Ru-Rl+1
        histogram(Sbs[j],bins=mbn,label="",title="Strains with $(j+Rl-1) reactions")
        plot!(xlabel="Metabolite used by reaction")
        savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/WhichSubs$(Rl)-$(Ru)$(syn)$(Ni)R=$(j+Rl-1).png")
        SbsT = cat(SbsT,Sbs[j],dims=1)
    end
    histogram(SbsT,bins=mbn,label="",title="All strains",xlabel="Metabolite used by reaction")
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/WhichSubs$(Rl)-$(Ru)$(syn)$(Ni)All.png")
    # Make plot for strength of generalism
    plot(title="Generalism vs Substrate Diversity ($(Rl)-$(Ru) $(syn))",ylabel="Average number of reactions")
    plot!(xlabel="Prob of viable reaction")
    for i = 1:M-1
        if fnd[i] == true && length(avR[i]) > 1
            scatter!([i/(M-1)],[mean(avR[i])],yerror=[std(avR[i])],label="")
        elseif fnd[i] == true
            scatter!([i/(M-1)],[mean(avR[i])],label="")
        end
    end
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/GenvsDiv$(Rl)-$(Ru)$(syn)$(Ni).png")
    return(nothing)
end

# function to take in multiple parameter sets and plot composite graphs
function multi_sets()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 3
        error("Need to provide number of parameter sets, number of repeats and initial number of strains")
    end
    # Preallocate the variables I want to extract from the input
    Ns = 0
    Nr = 0
    Ni = 0
    # Check that all arguments can be converted to integers
    try
        Ns = parse(Int64,ARGS[1])
        Nr = parse(Int64,ARGS[2])
        Ni = parse(Int64,ARGS[3])
    catch e
        error("number of parameter sets and repeats must both be integer")
    end
    # Preallocate vectors to store input
    Rls = zeros(Int64,Ns)
    Rus = zeros(Int64,Ns)
    syns = fill(false,Ns)
    ens = fill("",Ns)
    for i = 1:Ns
        vld = false
        while vld == false
            # Get input from the user
            println("Provide reaction lower bound")
            Rl = readline()
            try
                Rl = parse(Int64,Rl)
            catch e
                Rl = 0 # Set to invalid value
            end
            println("Provide reaction upper bound")
            Ru = readline()
            try
                Ru = parse(Int64,Ru)
            catch e
                Ru = 0 # Set to invalid value
            end
            println("Syntrophy true/false?")
            syn = readline()
            try
                syn = parse(Bool,syn)
            catch e
                syn = undef # Set to invalid value
            end
            println("Provide energy supply level (l/i/h)")
            en = readline()
            # Repeat checking step
            rpt = false
            # Loop over previous entries
            for j = 1:(i-1)
                if Ru == Rus[j] && Rl == Rls[j] && syn == syns[j] && en == ens[j]
                    rpt = true
                end
            end
            # Check if there is a directory containing these parameter sets
            if rpt == true
                println("Parameter set repeated, please enter an orginal one")
            elseif ~isdir("Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)")
                println("Parameter set doesn't exist please reenter")
            else
                vld = true
                Rls[i] = Rl
                Rus[i] = Ru
                syns[i] = syn
                ens[i] = en
            end
        end
    end
    # Container to store number of survivors
    svs = zeros(Int64,Ns,Nr)
    # Container to store metabolite diversity
    mbs = zeros(Int64,Ns,Nr)
    # Container to store lowest metabolite
    hmb = zeros(Int64,Ns,Nr)
    # Container to store mean abundances
    mna = zeros(Float64,Ns,Nr)
    # Container to store median abundances
    mda = zeros(Float64,Ns,Nr)
    # Container to store total abundances
    tab = zeros(Float64,Ns,Nr)
    # Empty container for abundances
    abs = Array{Array{Float64,1},2}(undef,Ns,Nr)
    # Preallocate labels
    lbs = Array{String}(undef,Ns)
    # Loop over parameter sets
    for i = 1:Ns
        for j = 1:Nr
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
            ps = load(pfile,"ps")
            inf_out = load(ofile,"inf_out")
            # Save number of survivors
            svs[i,j] = ps.N
            # Save abundances
            abs[i,j] = inf_out[1:ps.N]
            # Loop over metabolites to find those with non-zero concentrations
            cm = 0 # Set up counter
            mm = 0 # Lowest metabolite
            for k = 2:ps.M
                if inf_out[ps.N+k] > 0.0
                    # Increment counter
                    cm += 1
                    # And save new minimum metabolite
                    mm = k
                end
            end
            # Save results to vector
            mbs[i,j] = cm
            hmb[i,j] = mm
            # Save mean and median abundances
            mna[i,j] = mean(inf_out[1:ps.N])
            mda[i,j] = median(inf_out[1:ps.N])
            # Also save total abundance
            tab[i,j] = sum(inf_out[1:ps.N])
        end
        # Make and save label
        if Rls[i] != Rus[i]
            lbs[i] = "$(Rls[i])-$(Rus[i]) $(syns[i]) $(ens[i])"
        else
            lbs[i] = "$(Rus[i]) $(syns[i]) $(ens[i])"
        end
    end
    # Make empty containers to store
    tl = Array{String,1}(undef,Ns*Nr)
    tsv = Array{Int64,1}(undef,Ns*Nr)
    tmn = Array{Float64,1}(undef,Ns*Nr)
    tmd = Array{Float64,1}(undef,Ns*Nr)
    tdv = Array{Int64,1}(undef,Ns*Nr)
    tta = Array{Float64,1}(undef,Ns*Nr)
    enl = Array{String,1}(undef,Ns*Nr)
    # Fill out with data
    for i = 1:Ns
        for j = 1:Nr
            tl[(i-1)*Nr+j] = lbs[i]
            tsv[(i-1)*Nr+j] = svs[i,j]
            tmn[(i-1)*Nr+j] = mna[i,j]
            tmd[(i-1)*Nr+j] = mda[i,j]
            tdv[(i-1)*Nr+j] = mbs[i,j]
            tta[(i-1)*Nr+j] = tab[i,j]
            enl[(i-1)*Nr+j] = ens[i]
        end
    end
    # Collect everything into one data frame
    survivors = DataFrame(PSet=tl,ns=tsv,mn=tmn,md=tmd,sdv=tdv,ta=tta,en=enl)
    # Need to make a second data frame
    absT = Float64[]
    lbT = String[]
    enT = String[]
    mn_abs = zeros(Ns)
    sd_abs = zeros(Ns)
    # Fill out with data
    for i = 1:Ns
        # Make temporary vector for parameter set
        abst = Float64[]
        for j = 1:Nr
            # Cat data into a single vector
            abst = cat(abst,abs[i,j],dims=1)
            # Temporary vector of labels
            tl = fill(lbs[i],length(abs[i,j]))
            # Cat into long vector of labels
            lbT = cat(lbT,tl,dims=1)
            # Temporary vector of labels
            el = fill(ens[i],length(abs[i,j]))
            # Cat into long vector of energies
            enT = cat(enT,el,dims=1)
        end
        # Calculate and store mean of temporary vector
        mn_abs[i] = mean(abst)
        sd_abs[i] = std(abst)
        # Cat temporary vector into overall vector
        absT = cat(absT,abst,dims=1)
    end
    # Collect everything into one data frame
    abundances = DataFrame(PSet=lbT,abun=absT,en=enT)
    # Find indices of labels
    lbI = sortperm(lbs)
    # Find corresponding indices of reordered labels
    pos = zeros(Float64,Ns)
    # Find posistions for mean and std
    for i = 1:Ns
        p = 0.5
        if syns[i] == true
            p += 1.0
        end
        if ens[i] == "i"
            p += 2.4
        elseif ens[i] == "l"
            p += 4.8
        end
        pos[i] = p
    end
    # Setup plotting
    pyplot()
    theme(:wong2,dpi=200)
    wongc = get_color_palette(wong_palette,57)
    # Want to do the plotting here
    plot(title="Ecosystem diversity",ylabel="Number of surviving strains")
    @df survivors violin!(:PSet,:ns,linewidth=0,label="",color=wongc[2],group=:en)
    @df survivors boxplot!(:PSet,:ns,color=wongc[4],fillalpha=0.75,linewidth=2,label="",group=:en)
    # Plot means
    for i = 1:Ns
        # Calculate mean
        mn = mean(svs[i,:])
        sdn = std(svs[i,:])
        scatter!([pos[i]],[mn],yerror=[sdn],label="",shape=:star5,color=wongc[5],ms=10,msc=wongc[5])
    end
    savefig("Output/Diversity.png")
    plot(title="Mean abundances",ylabel="Mean abundance",yaxis=:log10)
    @df survivors violin!(:PSet,:mn,linewidth=0,label="",color=wongc[2],group=:en)
    @df survivors boxplot!(:PSet,:mn,color=wongc[4],fillalpha=0.75,linewidth=2,label="",group=:en)
    # Plot means
    for i = 1:Ns
        # Calculate mean
        mn = mean(mna[i,:])
        sdn = std(mna[i,:])
        scatter!([pos[i]],[mn],yerror=[sdn],label="",shape=:star5,color=wongc[5],ms=10,msc=wongc[5])
    end
    savefig("Output/MeanAbund.png")
    plot(title="Median abundances",ylabel="Median abundance",yaxis=:log10)
    @df survivors violin!(:PSet,:md,linewidth=0,label="",color=wongc[2],group=:en)
    @df survivors boxplot!(:PSet,:md,color=wongc[4],fillalpha=0.75,linewidth=2,label="",group=:en)
    # Plot means
    for i = 1:Ns
        # Calculate mean
        mn = mean(mda[i,:])
        sdn = std(mda[i,:])
        scatter!([pos[i]],[mn],yerror=[sdn],label="",shape=:star5,color=wongc[5],ms=10,msc=wongc[5])
    end
    savefig("Output/MedianAbund.png")
    plot(title="Substrate diversification",ylabel="Number of substrates")
    @df survivors violin!(:PSet,:sdv,linewidth=0,label="",color=wongc[2],group=:en)
    @df survivors boxplot!(:PSet,:sdv,color=wongc[4],fillalpha=0.75,linewidth=2,label="",group=:en)
    # Plot means
    for i = 1:Ns
        # Calculate mean
        mn = mean(mbs[i,:])
        sdn = std(mbs[i,:])
        scatter!([pos[i]],[mn],yerror=[sdn],label="",shape=:star5,color=wongc[5],ms=10,msc=wongc[5])
    end
    savefig("Output/SubDiv.png")
    plot(title="All abundances",ylabel="Strain abundance",yaxis=:log10)
    @df abundances violin!(:PSet,:abun,linewidth=0,label="",color=wongc[2],group=:en)
    @df abundances boxplot!(:PSet,:abun,color=wongc[4],fillalpha=0.75,linewidth=2,label="",group=:en)
    # Plot means
    for i = 1:Ns
        scatter!([pos[i]],[mn_abs[i]],yerror=[sd_abs[i]],label="",shape=:star5,color=wongc[5],ms=10,msc=wongc[5])
    end
    savefig("Output/AllAbund.png")
    plot(title="Total abundances",ylabel="Total abundance (per ecosystem)",yaxis=:log10)
    @df survivors violin!(:PSet,:ta,linewidth=0,label="",color=wongc[2],group=:en)
    @df survivors boxplot!(:PSet,:ta,color=wongc[4],fillalpha=0.75,linewidth=2,label="",group=:en)
    # Plot means
    for i = 1:Ns
        # Calculate mean
        mn = mean(tab[i,:])
        sdn = std(tab[i,:])
        scatter!([pos[i]],[mn],yerror=[sdn],label="",shape=:star5,color=wongc[5],ms=10,msc=wongc[5])
    end
    savefig("Output/TotalAbund.png")
    return(nothing)
end

# Function to plot survivorship with time
function plot_survivors()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 7
        error("Insufficent inputs provided (looking for 7)")
    end
    # Preallocate the variables I want to extract from the input
    Rl = 0
    Ru = 0
    syn = true
    rps = 0
    Ni = 0
    fT = 0.0
    sT = 0.0
    # Check that all arguments can be converted to integers
    try
        Rl = parse(Int64,ARGS[1])
        Ru = parse(Int64,ARGS[2])
        syn = parse(Bool,ARGS[3])
        rps = parse(Int64,ARGS[4])
        Ni = parse(Int64,ARGS[5])
        fT = parse(Float64,ARGS[6])
        sT = parse(Float64,ARGS[7])
    catch e
        error("need to provide 4 integers, a bool, and two floats")
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
        error("number of initial strains can't be less than 1")
    end
    if fT <= 0.0 || fT > 1.0
        error("final time can't be greater than 1, or less than zero")
    end
    if sT >= fT || sT < 0.0
        error("starting time can't be equal to final time, or less than zero")
    end
    println("Compiled!")
    # Want to look at all three energy conditions
    ens = ["l","i","h"]
    # Choosing to sample a thousand points for now
    ips = 1000
    # Read in first data file (chose l)
    ofile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)l/OutputReacs$(Rl)-$(Ru)Syn$(syn)Run1Ns$(Ni).jld"
    if ~isfile(ofile)
        error("run 1 energy supply l is missing an output file")
    end
    pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)l/ParasReacs$(Rl)-$(Ru)Syn$(syn)Run1Ns$(Ni).jld"
    if ~isfile(ofile)
        error("run 1 energy supply l is missing an parameter file")
    end
    ps = load(pfile,"ps")
    T = load(ofile,"T")
    # And extract maximum time
    Tmax = T[end]
    # Make vector of times to check at
    Ts = collect(range(Tmax*sT,Tmax*fT,length=ips))
    # Preallocate number of survivors, substrate diversities, total abundance
    svs = zeros(ips,rps,length(ens))
    dv = zeros(ips,rps,length(ens))
    tab = zeros(ips,rps,length(ens))
    # Set threshold for substrate being properly diversified
    tsh = 1e-7
    # Loop over repeats
    for i = 1:rps
        for j = 1:length(ens)
            # Read in relevant files
            ofile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(ens[j])/OutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ni).jld"
            if ~isfile(ofile)
                error("run $(i) energy supply $(ens[j]) is missing an output file")
            end
            # Load full dynamics
            C = load(ofile,"C")
            T = load(ofile,"T")
            # Loop over the time points
            for k = 1:ips
                # Find first time point greater than or equal to one were looking for
                ind = findfirst(x->x>=Ts[k],T)
                # If time points are equal just save number of survivors
                if T[ind] == Ts[k]
                    svs[k,i,j] = count(x->x>0.0,C[ind,1:Ni])
                    dv[k,i,j] = count(x->x>0.0,C[ind,(Ni+1):(Ni+ps.M-1)])
                    tab[k,i,j] = sum(C[ind,1:Ni])
                else
                    # Otherwise need to (linearly) interpolate
                    dT = (T[ind]-Ts[k])/(T[ind]-T[ind-1])
                    svs[k,i,j] = (1-dT)*count(x->x>0.0,C[ind,1:Ni]) + dT*count(x->x>0.0,C[ind-1,1:Ni])
                    dv[k,i,j] = (1-dT)*count(x->x>0.0,C[ind,(Ni+1):(Ni+ps.M-1)]) + dT*count(x->x>0.0,C[ind-1,(Ni+1):(Ni+ps.M-1)])
                    tab[k,i,j] = (1-dT)*sum(C[ind,1:Ni]) + dT*sum(C[ind,1:Ni])
                end
            end
        end
    end
    # Preallocate means and sds
    msvs = zeros(ips,length(ens))
    sdsvs = zeros(ips,length(ens))
    mdv = zeros(ips,length(ens))
    sddv = zeros(ips,length(ens))
    mta = zeros(ips,length(ens))
    sdta = zeros(ips,length(ens))
    # Find mean and standard errors of survivor numbers and substrate diversification
    for i = 1:ips
        for j = 1:length(ens)
            msvs[i,j] = mean(svs[i,:,j])
            sdsvs[i,j] = sem(svs[i,:,j])
            mdv[i,j] = mean(dv[i,:,j])
            sddv[i,j] = sem(dv[i,:,j])
            mta[i,j] = mean(tab[i,:,j])
            sdta[i,j] = sem(tab[i,:,j])
        end
    end
    # Set up plotting
    pyplot()
    # Set a color-blind friendly palette
    theme(:wong2,dpi=200)
    # Make labels
    lb = ["low","intermediate","high"]
    # Plot survivors
    plot(title="Survival ($(Rl)-$(Ru) $(syn))",xlabel="Time",ylabel="Number of surviving strains")
    for i = 1:length(ens)
        plot!(Ts,msvs[:,i],ribbon=sdsvs[:,i],label=lb[i])
    end
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)/SvTime$(Rl)-$(Ru)$(syn)$(Ni).png")
    # Plot substrate diversity
    plot(title="Diversification ($(Rl)-$(Ru) $(syn))",xlabel="Time",ylabel="Number of substrates")
    for i = 1:length(ens)
        plot!(Ts,mdv[:,i],ribbon=sddv[:,i],label=lb[i])
    end
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)/DvTime$(Rl)-$(Ru)$(syn)$(Ni).png")
    # Plot total abundances
    plot(title="Ecosystem abundance ($(Rl)-$(Ru) $(syn))",xlabel="Time",ylabel="Total abundance")
    for i = 1:length(ens)
        plot!(Ts,mta[:,i],ribbon=sdta[:,i],label=lb[i])
    end
    savefig("Output/$(Rl)-$(Ru)$(syn)$(Ni)/TotalAbTime$(Rl)-$(Ru)$(syn)$(Ni).png")
    return(nothing)
end

if length(ARGS) == 3
    @time multi_sets()
elseif length(ARGS) == 7
    @time plot_survivors()
else
    @time basic_info()
end
