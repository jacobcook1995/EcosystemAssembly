# Script to test outcome of just the number of ribosomes trade-off
using TradeOff
using Random
using LaTeXStrings
using Plots
import PyPlot

# Function to make the kind of microbe needed to test the ω tradeoff
function new_mic_ωt(M::Int64,ω::Float64,KΩ::Float64,Kγ::Float64,Rl::Int64,Ru::Int64,ps::TOParameters,d::Float64)
    # First generate random unique indetifier for this pool
    PID = randstring(['0':'9'; 'a':'f'])
    # Only generating one microbe so will be ID: 1
    ID = 1
    # Print out that this is happening
    println("Generating random microbe with identifer: $(PID)")
    # Assume that half saturation occurs at a quarter κ/δ
    KS = (1/4)*5.5e-3
    # Arbitary number that seems to give decent survival
    kc = 10.0
    # The reversibility factor remains the same as previously
    kr = 10.0
    # Use formula to calculate how many reactions are implied
    O = 2*M - 3
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # Cell mass is taken from Bremer H, Dennis P (1996) Modulation of chemical
    # composition and other parameters of the cell by growth rate (Book chapter).
    MC = 10^8
    # Elongation rate also taken from Bremer (1996), convert from minutes to seconds
    γm = 1260.0/60.0
    # Now preallocate protein masses
    n = zeros(Int64,3)
    # ribosome mass taken from Keseler IM, et al. (2011)
    n[1] = 7459
    # Other protein mass averaged from Brandt F, et al. (2009)
    n[2:3] .= 300
    # The number of ATP per translation step, including the cost of amino acid sythesis
    # This figure is taken from Lynch and Marinov 2015
    ρ = 27.55
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # Number of doublings required to dilute to 1%
    fd = log(100)/log(2)
    # For each microbe generate random set of reactions
    R, Reacs = choose_reactions(O,Rl,Ru)
    # Make vectors of the (fixed) kinetic parameters
    kcs = kc.*ones(R)
    KSs = KS.*ones(R)
    krs = kr.*ones(R)
    # Reactions given random proportional weightings, done this in the simplest way possible
    ϕP = rand(R)
    ϕP = ϕP/sum(ϕP)
    # Find corresponding η's for these reactions
    η = zeros(R)
    # Set ratio for equilbrium
    mratio = 1e5
    # Loop over all reactions
    for i = 1:length(η)
        # Identify which reaction we are considering
        I = Reacs[i]
        # Find corresponding Gibbs free energy change
        dG = ps.reacs[I].ΔG0
        # And then use to determine an upper bound on η
        η[i] = -(dG + Rgas*T*log(mratio))/(ΔGATP)
    end
    # Can finally generate microbe
    mic = make_Microbe(MC,γm,ρ,Kγ,Pb,d,ϕH,KΩ,fd,ω,R,Reacs,η,kcs,KSs,krs,n,ϕP,ID,PID)
    return(mic)
end

# function to look at competition between two strains
function comp_2()
    # Two metabolites 1 reaction for the sake of simplicity
    M = 2
    O = 1
    # Somewhat arbitary free energy range
    μrange = 4e5
    # First generate parameter set
    ps = initialise(M,O,μrange)
    # Only one reaction so upper and lower bound have to be one
    Rl = 1
    Ru = 1
    # Set parameters that vary
    ωs = [0.75,1.00]
    KΩ = 1e9
    Kγ = 2.5e8
    # Also allowing dilution rate to vary
    d = 6.0e-5
    # Preallocate array of two microbes
    mics = Array{Microbe,1}(undef,2)
    for i = 1:length(mics)
        mics[i] = new_mic_ωt(M,ωs[i],KΩ,Kγ,Rl,Ru,ps,d)
    end
    # Set reasonable time window
    Tmax = 1e8
    # Choose initial condition
    pop = 1000.0
    conc = 1e-5
    as = 1e5
    ϕs = 0.128
    # Then simulate
    C, T = doub_pop(ps,pop,conc,as,ϕs,mics,Tmax)
    # Finally plot this
    pyplot(dpi=200)
    # Considering 2 strains at the moment
    totN = 2
    # Plot all the populations
    p1 = plot(yaxis=:log10,ylabel="Population (# cells)",xlabel="Time")
    p2 = plot(ylabel="ATP",xlabel="Time")
    p3 = plot(ylabel=L"\phi_R",xlabel="Time")
    for i = 1:totN
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0)
        plot!(p1,T[inds],C[inds,i],label="ω = $(ωs[i])")
        plot!(p2,T,C[:,(totN+ps.M+i)],label="ω = $(ωs[i])")
        plot!(p3,T,C[:,(2*totN+ps.M+i)],label="ω = $(ωs[i])")
    end
    # plot horizontal line for half saturation
    hline!(p2,[Kγ],style=:dash,label="")
    # Then save figures
    savefig(p1,"Output/pops.png")
    savefig(p2,"Output/as.png")
    savefig(p3,"Output/fracs.png")
    # Plot all the concentrations
    p4 = plot(yaxis=:log10,ylabel="Concentration",xlabel="Time")
    for i = 1:ps.M
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,totN+i] .> 0)
        plot!(p4,T[inds],C[inds,totN+i],label="Metabolite $i")
    end
    savefig(p4,"Output/concs.png")
    return(nothing)
end

@time comp_2()
