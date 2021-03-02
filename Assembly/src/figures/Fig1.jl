# Script to plot elements needed for figure 1
using Assembly
using Plots
using LaTeXStrings
using JLD
import PyPlot

# function to plot population dynamics
function popdyn(Rl::Int64,Ru::Int64,syn::Bool,Nr::Int64,Ns::Int64,en::String)
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
    # Find maximum time
    Tmax = T[end]
    # Now move onto plotting
    pyplot()
    theme(:wong2,dpi=200)
    # Plot all the populations
    plot(title="Population dynamics",yaxis=:log10,xlabel="Time (s)",ylabel="Population (# cells)")
    for i = 1:Ns
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0) .& (T .< Tmax/4)
        plot!(T[inds],C[inds,i],label="")
    end
    savefig("Output/Fig1/fullpops.png")
    return(nothing)
end

# Function to make plots to show the shape of the thermodynamic tradeoff
function plt_trdff(Rl::Int64,Ru::Int64,syn::Bool,runN::Int64,en::String,Ni::Int64)
    println("Compiled!")
    # Read in relevant files
    pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(runN)Ns$(Ni).jld"
    if ~isfile(pfile)
        error("run $(runN) is missing a parameter file")
    end
    ofile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/RedOutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(runN)Ns$(Ni).jld"
    if ~isfile(ofile)
        error("run $(runN) is missing an output file")
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
    # Loop over η values
    for i = 1:length(ηs)
        # Calculate thermodynamic inhibition
        θs[i] = θ_smooth(S,P,ps.T,ηs[i],ΔG)
        # Then use to calculate rate
        rs[i] = qs(ps.mics[1],S,P,E,θs[i])
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
    savefig("Output/Fig1/TrdOff.png")
    # Do plot of just the tradeoff
    plot(ηs,as,label="",xlabel="ATP per reaction event",ylabel="ATP production rate")
    savefig("Output/Fig1/BareTrdOff.png")
    # Want a detailed visulisation of the peak
    inds = findall(x->(5.25<=x<=5.75),ηs)
    plot(ηs[inds],rs[inds],label="Reaction rate")
    plot!(ηs[inds],θs[inds],label="Inhibition")
    plot!(ηs[inds],as2[inds],label="ATP rate",xlabel="ATP per reaction event")
    savefig("Output/Fig1/PeakTrdOff.png")
    return(nothing)
end

# @time popdyn(1,5,true,78,250,"l")
@time plt_trdff(1,5,true,67,"i",250)
