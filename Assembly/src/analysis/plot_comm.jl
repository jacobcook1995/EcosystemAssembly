# This script exists to plot the output of form_comm.jl
using Assembly
using Plots
using LaTeXStrings
using JLD
using SymPy
import PyPlot

# function to calculate the dissipation for an assembled ecosystem
function dissipation(ps::FullParameters,out::Array{Float64,1})
    # check that parameter set is sensible given the output
    if length(out) != ps.M + 3*ps.N
        error("parameter set doesn't match output")
    end
    # Set dissipation to zero
    dsp = 0
    # Loop over number of strains
    for i = 1:ps.N
        # Isolate this strain
        mic = ps.mics[i]
        # Loop over reactions of this strain
        for j = 1:mic.R
            # Find appropriate reaction
            r = ps.reacs[mic.Reacs[j]]
            # Find amount of energy that this reaction dissipates
            Fd = -(r.ΔG0 + Rgas*ps.T*log(out[ps.N+r.Prd]/out[ps.N+r.Rct]) + mic.η[j]*ΔGATP)
            # Find amount of enzyme E
            E = Eα(out[2*ps.N+ps.M+i],mic,j)
            # Then find the rate that this reaction proceeds at
            q = qs(out[ps.N+r.Rct],out[ps.N+r.Prd],E,j,mic,ps.T,r)
            # Check if reaction actually occurs
            if q != 0.0
                dsp += q*Fd
            end
        end
    end
    return(dsp)
end

# A function to read in saved data and interpret it
function interpret()
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
    # Preallocate memory to store the number of survivors in
    srv = zeros(rps)
    # Preallocate vector to store percentage of free energy dissipated (under standard conditions)
    pd = []
    # Preallocate vector of dissipation rates
    dsp = zeros(rps)
    # Set count of stable systems to zero
    cnt = 0
    # Preallocate vector to store stability times
    sT = []
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
        if ~isfile(ofile)
            error("run $(i) is missing an extinct file")
        end
        # Then load in data
        ps = load(pfile,"ps")
        out = load(ofile,"out")
        # Save number of surviving species
        srv[i] = ps.N
        # Store all percentage dissipations
        for j = 1:ps.N
            # Find vector of η values
            ηs = ps.mics[j].η
            # make temporary vector for percentages
            pdt = zeros(length(ηs))
            # loop over η values
            for k = 1:length(ηs)
                # Find standard Gibbs free energy
                dG = ps.reacs[ps.mics[j].Reacs[k]].ΔG0
                # Calculate percentage dissipated (under standard conditions)
                pdt[k] = -dG/(ηs[k]*ΔGATP)
            end
            pd = cat(pd,pdt,dims=1)
        end
        # Remove any negative metabolite concentrations
        for i = (ps.N+1):(ps.N+ps.M)
            if out[i] < 0.0
                out[i] = 0.0
            end
        end
        # Now its safe to calculate the dissipation
        dsp[i] = dissipation(ps,out)
    end
    # Now move onto plotting
    println("Data read in")
    pyplot(dpi=200)
    # Save number of surviving species as a histogram
    histogram(srv,label="")
    plot!(title="$(R) reactions",xlabel="Number of surving strains")
    savefig("Output/Type$(R)SvHist.png")
    histogram(pd,label="")
    plot!(title="$(R) reactions",xlabel="Percentage of free energy dissipated")
    savefig("Output/Type$(R)ηHist.png")
    histogram(dsp,label="")
    plot!(title="$(R) reactions",xlabel="Dissipation rate")
    savefig("Output/Type$(R)DispHist.png")
    # Read in stability time data
    sT = load(tfile,"sT")
    Ns = length(sT) - sum(isnan.(sT))
    # Then plot histogram
    histogram(sT,label="$(Ns)/$(length(sT)) stable")
    plot!(title="$(R) reactions",xlabel="Time for system to stabilize")
    savefig("Output/Type$(R)StabHist.png")
    # NEXT I WANT TO PLOT ILLUSTRATIVE EXAMPLES OF THE 
    return(nothing)
end

@time interpret()

# SAVE THIS FOR LATER
# # Useful to also have the extinction data available
# ded = load(efile,"ded")
# # Setup population plot
# pyplot(dpi=200)
# p1 = plot(xlabel="Time",ylabel="Population",yaxis=:log10)
# # Find number of ignored microbes
# N = ps.N + length(ded)
# for i = 1:N
#     # Find and eliminate zeros so that they can be plotted on a log plot
#     inds = (C[:,i] .> 0)
#     plot!(p1,T[inds],C[inds,i],label="")
# end
# savefig(p1,"Output/PopvsTime.png")
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
