# This script exists to plot the output of form_comm.jl
using Assembly
using Plots
using LaTeXStrings
using JLD
using SymPy
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
    # Check that stability data exists
    tfile = "Data/Type$(R)/StabTimesType$(R).jld"
    if ~isfile(tfile)
        error("no stability time file")
    end
    println("Compiled!")
    # Now move onto plotting
    pyplot(dpi=200)
    # Read in stability time data
    sT = load(tfile,"sT")
    Ns = length(sT) - sum(isnan.(sT))
    # Then plot histogram
    histogram(sT,label="$(Ns)/$(length(sT)) stable")
    plot!(title="$(R) reactions",xlabel="Time for system to stabilize")
    savefig("Output/Type$(R)StabHist.png")
    # Check if there are non-stable cases and plot one
    if Ns != length(sT)
        # Pick random number to select parameter set
        r = rand(1:(length(sT)-Ns))
        ns = findall(isnan,sT)
        # Read in relevant files
        pfile = "Data/Type$(R)/ParasType$(R)Run$(ns[r]).jld"
        if ~isfile(pfile)
            error("run $(ns[r]) is missing a parameter file")
        end
        ofile = "Data/Type$(R)/OutputType$(R)Run$(ns[r]).jld"
        if ~isfile(ofile)
            error("run $(ns[r]) is missing an output file")
        end
        efile = "Data/Type$(R)/ExtinctType$(R)Run$(ns[r]).jld"
        if ~isfile(efile)
            error("run $(ns[r]) is missing an extinct file")
        end
        ps = load(pfile,"ps")
        C = load(ofile,"C")
        T = load(ofile,"T")
        ded = load(efile,"ded")
        # Setup population plot
        p1 = plot(xlabel="Time",ylabel="Population",yaxis=:log10)
        # Find number of ignored microbes
        N = ps.N + length(ded)
        for i = 1:N
            # Find and eliminate zeros so that they can be plotted on a log plot
            inds = (C[:,i] .> 0)
            plot!(p1,T[inds],C[inds,i],label="")
        end
        savefig(p1,"Output/UnstabType$(R).png")
    end
    # Check if there are stable cases and plot one
    if Ns != 0
        # Pick random number to select parameter set
        r = rand(1:Ns)
        ns = findall(!isnan,sT)
        # Read in relevant files
        pfile = "Data/Type$(R)/ParasType$(R)Run$(ns[r]).jld"
        if ~isfile(pfile)
            error("run $(ns[r]) is missing a parameter file")
        end
        ofile = "Data/Type$(R)/OutputType$(R)Run$(ns[r]).jld"
        if ~isfile(ofile)
            error("run $(ns[r]) is missing an output file")
        end
        efile = "Data/Type$(R)/ExtinctType$(R)Run$(ns[r]).jld"
        if ~isfile(efile)
            error("run $(ns[r]) is missing an extinct file")
        end
        ps = load(pfile,"ps")
        C = load(ofile,"C")
        T = load(ofile,"T")
        ded = load(efile,"ded")
        # Setup population plot
        p1 = plot(xlabel="Time",ylabel="Population",yaxis=:log10)
        # Find number of ignored microbes
        N = ps.N + length(ded)
        for i = 1:N
            # Find and eliminate zeros so that they can be plotted on a log plot
            inds = (C[:,i] .> 0)
            plot!(p1,T[inds],C[inds,i],label="")
        end
        vline!([sT[ns[r]]],label="",color=:red)
        savefig(p1,"Output/StabPopType$(R).png")
    end
    # Mark on these the time point of stabilization
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
    # Preallocate vector to store percentage of free energy dissipated (under standard conditions)
    pd = fill(Float64[],length(Ts))
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
            # Remove any negative metabolite concentrations
            for i = 1:ps.M
                if concs[i] < 0.0
                    concs[i] = 0.0
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
                    # Find vector of η values
                    ηs = ms[k].η
                    # make temporary vector for percentages
                    pdt = zeros(length(ηs))
                    # loop over η values
                    for l = 1:length(ηs)
                        # Find standard Gibbs free energy
                        dG = ps.reacs[ms[k].Reacs[l]].ΔG0
                        # Calculate percentage dissipated (under standard conditions)
                        pdt[l] = -dG/(ηs[l]*ΔGATP)
                    end
                    # Here I cat into the preallocate array of arrays
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
    pyplot(dpi=200)
    # Preallocate vector of labels
    Tlbs = Array{String,2}(undef,1,length(Ts))
    # Insert all elements to it
    for i = 1:length(Ts)
        Tlbs[i] = "T = $(Ts[i])s"
    end
    # Now plot the histograms
    m1 = L"^{-1}"
    histogram(srv,labels=Tlbs,fillalpha=0.75)
    plot!(title="$(R) reactions",xlabel="Number of surviving strains")
    savefig("Output/SvType$(R).png")
    histogram(dsp,labels=Tlbs,fillalpha=0.75)
    plot!(title="$(R) reactions",xlabel="Dissipation rate (Js$(m1))")
    savefig("Output/DispType$(R).png")
    histogram(b_dsp,labels=Tlbs,fillalpha=0.75)
    plot!(title="$(R) reactions",xlabel="Dissipation rate (Js$(m1)) per cell")
    savefig("Output/BDispType$(R).png")
    histogram(pd,labels=Tlbs,fillalpha=0.75)
    plot!(title="$(R) reactions",xlabel="Percentage of free energy dissipated")
    savefig("Output/RDispType$(R).png")
    Find and plot percentage active
    histogram(100*(1 .- inact./(inact.+ob.+nob)),labels=Tlbs,fillalpha=0.75)
    plot!(title="$(R) reactions",xlabel="Percentage of reactions active")
    savefig("Output/PerActType$(R).png")
    histogram(100*ob./(ob.+nob),labels=Tlbs,fillalpha=0.75)
    plot!(title="$(R) reactions",xlabel="Percentage of reactions obligate")
    savefig("Output/PerObType$(R).png")
    return(nothing)
end

# At the moment want to run both plotting scripts in succession
@time stab_plots()
@time hist_time()

# SAVE THIS FOR LATER
# # Useful to also have the extinction data available
# ded = load(efile,"ded")

# # Now plot the metabolites
# plot(T,C[:,(N+1):(N+ps.M)],xlabel="Time",label="",ylabel="Concentration")
# savefig("Output/MetabolitevsTime.png")
# # Now plot the energy concentrations
# plot(T,C[:,(N+ps.M+1):(2*N+ps.M)],xlabel="Time",label="",ylabel="Cell energy conc")
# savefig("Output/EnergyvsTime.png")
# plot(T,C[:,2*N+ps.M+1:end],xlabel="Time",label="",ylabel=L"\phi_R")
# savefig("Output/FractionvsTime.png")
# # Calculate reaction inhibition
# # Find surviving microbes
# ms = []
# for i = 1:N
#     if C[end,i] > 1e-8
#         ms = cat(ms,i,dims=1)
#     end
# end
# # Make array of bools to store results
# # To generalise this I would have to figure out how to make
# # an array with rows of differing lengths
# rs = fill(false,(length(ms),ps.mics[1].R))
# # Then for these microbes find viable reactions
# for i = 1:length(ms)
#     # Find corresponding microbes
#     mic = ps.mics[ms[i]]
#     for j = 1:mic.R
#         # find reaction
#         r = ps.reacs[mic.Reacs[j]]
#         if C[end,N+r.Rct] > 1e-8 && C[end,N+r.Prd] > 1e-8
#             rs[i,j] = true
#         end
#     end
# end
# # Preallocate array of θs, this also will take some work to generalise
# θs = zeros(length(C[:,1]),size(rs,1)*size(rs,2))
# for i = 1:length(ms)
#     for j = 1:ps.mics[1].R
#         # Find appropriate reaction
#         r = ps.reacs[ps.mics[ms[i]].Reacs[j]]
#         for k = 1:length(C[:,1])
#             θs[k,(i-1)*ps.mics[1].R+j] = θ_smooth(C[k,N+r.Rct],C[k,N+r.Prd],ps.T,ps.mics[ms[i]].η[j],r.ΔG0)
#         end
#     end
# end
# # Then plot inhibitions
# plot(T,θs,xlabel="Time",label="",ylabel=L"\theta")
# savefig("Output/InhibvsTime.png")
