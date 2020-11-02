# This script exists to plot the output of form_comm.jl
using Assembly
using StatsPlots
using LaTeXStrings
using JLD
using StatsBase
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
function f_ext_un(x::StepRangeLen)
    # Preallocate output
    y = zeros(length(x))
    for i = 1:length(x)
        # Check if it's high enough to be in first uniform range
        if x[i] <= 0.8177
            y[i] += 1
        end
        # Check if it's high enough to be in second uniform range
        if x[i] <= 0.9089
            y[i] += 1
        end
        # Then subtract if too low
        if x[i] < -0.0409
            y[i] -= 1
        end
        if x[i] < -0.0818
            y[i] -= 1
        end
    end
    return(y)
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
    # To store number of ecosystems that reaction is present in
    pres = zeros(ps.O)
    # Data structure to store abundances
    abnd = zeros(ps.O,1)
    # Loop over repeats
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
            # Use to calculate percentage dissipated (under standard conditions)
            efT = (ps.mics[j].η*ΔGATP.+dG)./(dG)
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
            ef[j] = (ps.mics[j].η[ind]*ΔGATP+dG)/(dG)
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
    # Plot "presence" data
    groupedbar([nR.-pres pres],bar_position=:stack,labels=["Absent" "Present"])
    plot!(xticks=(1:ps.O,xs),xlabel="Reaction",ylabel="Number of ecosystems")
    plot!(title="$(R) reactions per strain",legend=:outerright)
    savefig("Output/Type$(R)/PresenseType$(R).png")
    # Now plot "abundance" data
    lbs = Array{String,2}(undef,1,size(abnd,2))
    for i = 1:length(lbs)
        lbs[i] = string(i-1)
    end
    groupedbar(abnd,bar_position=:stack,labels=lbs,legend=:outerright)
    plot!(xticks=(1:ps.O,xs),xlabel="Reaction",ylabel="Number of ecosystems")
    plot!(title="$(R) reactions per strain")
    savefig("Output/Type$(R)/AbundanceType$(R).png")
    # Plot histogram of the number of active reactions
    histogram(cA,bins=range(0,stop=O+1,length=O+2))
    plot!(title="$(R) reactions per strain",xlabel="Number of active reactions")
    savefig("Output/Type$(R)/ActiveReactionsType$(R).png")
    histogram(mf,bins=range(1,stop=O+1,length=O+1),xticks=(1:ps.O,xs))
    plot!(title="$(R) reactions per strain",xlabel="Reaction with greatest flux")
    savefig("Output/Type$(R)/MaxFluxType$(R).png")
    # Now find average flux
    fRT /= nR
    # Then plot
    bar(fRT,xticks=(1:ps.O,xs),label="",xlabel="Reaction",ylabel="Average flux")
    plot!(title="$(R) reactions per strain")
    savefig("Output/Type$(R)/AverageReacsType$(R).png")
    # Then find and plot average mass renormalised flux
    fRTm /= nR
    bar(fRTm,xticks=(1:ps.O,xs),label="",xlabel="Reaction",ylabel="Mass renormalised average flux")
    plot!(title="$(R) reactions per strain")
    savefig("Output/Type$(R)/AverageReacsMassType$(R).png")
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
    plot!(title="$(R) reactions per strain",xlabel="Reaction",xticks=(1:ps.O,xs))
    savefig("Output/Type$(R)/PercentFluxType$(R).png")
    bar(perc2,label="",ylabel="Percentange of flux (ignoring singles)")
    plot!(title="$(R) reactions per strain",xlabel="Reaction",xticks=(1:ps.O,xs))
    savefig("Output/Type$(R)/RedPercentFluxType$(R).png")
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
    bins5 = range(-0.1,stop=1.00,length=500)
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
    plot!(title="$(R) reactions per strain",xlabel=L"K_S")
    savefig("Output/Type$(R)/TopKSType$(R).png")
    bar(h2,label="Main reaction")
    μ2 = log(10.0)
    σ2 = log(2)
    d2 = f_ln(bins2,μ2,σ2)
    d2 /= maximum(d2)
    plot!(bins2,maximum(h2.weights)*d2,label="Orginal distribution")
    plot!(title="$(R) reactions per strain",xlabel=L"k_c")
    savefig("Output/Type$(R)/TopkcType$(R).png")
    bar(h3,label="Main reaction")
    μ3 = log(10.0)
    σ3 = log(2)
    d3 = f_ln(bins3,μ3,σ3)
    d3 /= maximum(d3)
    plot!(bins3,maximum(h3.weights)*d3,label="Orginal distribution")
    plot!(title="$(R) reactions per strain",xlabel=L"k_r")
    savefig("Output/Type$(R)/TopkrType$(R).png")
    bar(h4,label="Main reaction")
    μ4 = 1/R
    # No variation unless R > 1
    σ4 = 0.0
    if R != 1
        σ4 = (1/(R*sqrt(3)))*sqrt(1+(1/R))
    end
    d4 = f_nm(bins4,μ4,σ4)
    d4 /= maximum(d4)
    plot!(bins4,maximum(h4.weights)*d4,label="Orginal distribution")
    plot!(title="$(R) reactions per strain",xlabel=L"\phi_p")
    savefig("Output/Type$(R)/TopPhiPType$(R).png")
    bar(h5,label="Main reaction")
    d5 = f_ext_un(bins5)
    d5 /= maximum(d5)
    plot!(bins5,maximum(h5.weights)*d5,label="Orginal distribution",legend=:topright)
    plot!(title="$(R) reactions per strain",xlabel="Fraction of free energy dissipated")
    savefig("Output/Type$(R)/TopEfficencyType$(R).png")
    # Repeat process for full kinetics data
    bins1 = range(1e-4,stop=maximum(KSF),length=500)
    h1 = fit(Histogram,KSF,bins1,closed=:right)
    bins2 = range(1e-4,stop=maximum(kcF),length=500)
    h2 = fit(Histogram,kcF,bins2,closed=:right)
    bins3 = range(1e-4,stop=maximum(krF),length=500)
    h3 = fit(Histogram,krF,bins3,closed=:right)
    bins4 = range(1e-4,stop=maximum(ϕpF),length=500)
    h4 = fit(Histogram,ϕpF,bins4,closed=:right)
    bins5 = range(-0.1,stop=1.00,length=500)
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
    plot!(title="$(R) reactions per strain",xlabel=L"K_S")
    savefig("Output/Type$(R)/AllKSType$(R).png")
    bar(h2,label="All reactions")
    μ2 = log(10.0)
    σ2 = log(2)
    d2 = f_ln(bins2,μ2,σ2)
    d2 /= maximum(d2)
    plot!(bins2,maximum(h2.weights)*d2,label="Orginal distribution")
    plot!(title="$(R) reactions per strain",xlabel=L"k_c")
    savefig("Output/Type$(R)/AllkcType$(R).png")
    bar(h3,label="All reactions")
    μ3 = log(10.0)
    σ3 = log(2)
    d3 = f_ln(bins3,μ3,σ3)
    d3 /= maximum(d3)
    plot!(bins3,maximum(h3.weights)*d3,label="Orginal distribution")
    plot!(title="$(R) reactions per strain",xlabel=L"k_r")
    savefig("Output/Type$(R)/AllkrType$(R).png")
    bar(h4,label="All reactions")
    μ4 = 1/R
    # No variation unless R > 1
    σ4 = 0.0
    if R != 1
        σ4 = (1/(R*sqrt(3)))*sqrt(1+(1/R))
    end
    d4 = f_nm(bins4,μ4,σ4)
    d4 /= maximum(d4)
    plot!(bins4,maximum(h4.weights)*d4,label="Orginal distribution")
    plot!(title="$(R) reactions per strain",xlabel=L"\phi_p")
    savefig("Output/Type$(R)/AllPhiPType$(R).png")
    bar(h5,label="All reactions")
    d5 = f_ext_un(bins5)
    d5 /= maximum(d5)
    plot!(bins5,maximum(h5.weights)*d5,label="Orginal distribution",legend=:topright)
    plot!(title="$(R) reactions per strain",xlabel="Fraction of free energy dissipated")
    savefig("Output/Type$(R)/AllEfficencyType$(R).png")
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

# function I can use to test particular parameter sets
function test()
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
    # Read in relevant files
    pfile = "Data/Type$(R)/ParasType$(R)Run$(rps).jld"
    if ~isfile(pfile)
        error("run $(rps) is missing a parameter file")
    end
    ofile = "Data/Type$(R)/OutputType$(R)Run$(rps).jld"
    if ~isfile(ofile)
        error("run $(rps) is missing an output file")
    end
    efile = "Data/Type$(R)/ExtinctType$(R)Run$(rps).jld"
    if ~isfile(efile)
        error("run $(rps) is missing an extinct file")
    end
    ps = load(pfile,"ps")
    C = load(ofile,"C")
    T = load(ofile,"T")
    out = load(ofile,"out")
    ded = load(efile,"ded")
    # Find index I'm interested in
    ind = findfirst(x->x==out[4],C[end,:])
    m = ps.mics[4]
    # Now move onto plotting
    pyplot()
    theme(:wong2,dpi=200)
    plot(yaxis=:log10)
    # Find indicies
    is = zeros(Int64,ps.N)
    for i = 1:ps.N
        is[i] = findfirst(x->x==out[i],C[end,:])
    end
    N = ps.N + length(ded)
    # Something else
    for i = is
         # Find and eliminate zeros so that they can be plotted on a log plot
         inds = (C[:,i] .> 0)
         plot!(T[inds],C[inds,i],label="")
    end
    savefig("Output/testpop.png")
    plot(T,C[:,(N+1):(N+ps.M)])
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
    histogram(fl,label="",title="$(R) reactions per strain",xlabel="Flux")
    savefig("Output/Type$(R)/FluxesType$(R).png")
    histogram(flm,label="",title="$(R) reactions per strain",xlabel="Mass specific flux")
    savefig("Output/Type$(R)/FluxesMType$(R).png")
    histogram(Afl,label="",title="$(R) reactions per strain",xlabel="ATP flux")
    savefig("Output/Type$(R)/ATPFType$(R).png")
    histogram(Aflm,label="",title="$(R) reactions per strain",xlabel="Mass specific ATP flux")
    savefig("Output/Type$(R)/ATPFMType$(R).png")
    return(nothing)
end

@time test()
