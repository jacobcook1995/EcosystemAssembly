# Script to plot elements needed for figure 1
using Assembly
using Plots
using JLD
import PyPlot

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
    C = load(ofile,"C")
    T = load(ofile,"T")
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
    Ps = [10.0*P,P/10.0]
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
    wongc = wong2_palette()
    plot(ηs,rs,label="Reaction rate")
    plot!(ηs,θs,label="Inhibition")
    plot!(ηs,as2,label="ATP rate",xlabel="ATP per reaction event")
    savefig("Output/Fig1/TrdOff.png")
    # Do plot of just the tradeoff
    plot(ηs,as,label="",xlabel="ATP per reaction event",ylabel="ATP production rate")
    savefig("Output/Fig1/BareTrdOff.png")
    # Make a vector of η*rate
    as3 = ηs.*rs2
    # Define labels for the plot
    lbs = Array{String,2}(undef,1,2)
    lbs[1] = "High product concentration"
    lbs[2] = "Low product concentration"
    # Now calculate and plot syntrophy stuff
    plot(ηs,as3,xlabel="ATP per reaction event",ylabel="ATP production rate",labels=lbs,legend=:topleft)
    plot!(legendfontsize=12,guidefontsize=14,tickfontsize=10)
    # Add arrow between the two lines
    quiver!([5.55],[3e5],quiver=([-0.135],[0.0]),color=:red)
    savefig("Output/Fig1/SynTrdOff.png")
    # Find indicies of surviving strains
    is = zeros(Int64,ps.N)
    for i = 1:ps.N
        is[i] = findfirst(x->x==out[i],C[end,:])
    end
    # Only want to look at the very early time window
    Tend = 1e4
    # Define plot
    p1 = plot(ylabel="Population (# cells)",xlabel="Time (s)")
    for i = is
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (T .<= Tend)
        plot!(p1,T[inds],C[inds,i],label="")
    end
    savefig(p1,"Output/Fig1/EarlyDynamics.png")
    return(nothing)
end

@time plt_trdff(1,5,true,67,"i",250)
