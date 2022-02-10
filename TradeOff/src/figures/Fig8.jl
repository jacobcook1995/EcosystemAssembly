# Script to plot figure 8 (which shows parametrisation of Ω_1/2)
using TradeOff
using Plots
using Random
using ColorSchemes
using LaTeXStrings
using Plots.PlotMeasures
import PyPlot

# Function to plot the 7th figure (e.g. development without immigration)
function figure8()
    println("Compiled")
    # Define basic simulation parameters
    M = 2
    d = 6e-5
    μrange = 5e7 # Set high to avoid thermodynamic inhibition
    # Predefine a set of energy avaliabilities (η used as a proxy)
    ηs = collect(range(1.0, 20.0, length=20))
    # Predefine a set of omega values
    ωs = collect(range(0.1, 0.9, length=5))
    # Set a fixed value for Ω
    Ωf = 1e9
    # Define all the other fixed parameters
    PID = randstring(['0':'9'; 'a':'f'])
    MC = 10^8
    n = zeros(Int64,3)
    n[1] = 7459
    n[2:3] .= 300
    γm = 1260.0/60.0
    Kγ = 5e8
    χl = 29.0
    χu = 1e-20 # Turns off differences in efficency
    Pb = 0.7
    ϕH = 0.45
    fd = log(100)/log(2)
    KS = (1/4)*5.5e-3
    kc = 10.0
    kr = 10.0
    # Use formula to calculate how many reactions are implied
    O = 2*M - 3
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # Make parameter set
    ps = initialise(M,O,μrange)
    # Generate fixed reaction
    RP, ΔG = fix_reactions(O,M,μrange,T)
    # Preallocate vector to store single reaction
    reacs = Array{Reaction,1}(undef,O)
    for i = 1:O
        reacs[i] = make_Reaction(i,RP[i,1],RP[i,2],ΔG[i])
    end
    # Only one reaction, which all species possess
    R = 1
    Reacs = [1]
    # Preallocate array of fixed microbes
    fix = Array{Microbe,2}(undef,length(ηs),length(ωs))
    # Generate one species for each energy avaliabilities
    for i = 1:length(ηs)
        # And one for each ω parameter
        for j = 1:length(ωs)
            # Can finally generate microbe
            fix[i,j] = make_Microbe(MC,γm,Kγ,χl,χu,Pb,d,ϕH,Ωf,fd,ωs[j],R,Reacs,
                       [ηs[i]],[kc],[KS],[kr],n,[1.0],j+(i-1)*length(ηs),PID)
        end
    end
    # Preallocate vector of varying Ωs
    Ωvs = zeros(length(ωs))
    # Determine each Ω value based on corresponding ω
    for i = 1:length(ωs)
        Ωvs[i] = 1.5*(ωs[i]-0.1)*Ωf + 1*Ωf
    end
    # Preallocate array of varying microbes
    var = Array{Microbe,2}(undef,length(ηs),length(ωs))
    # Generate one species for each energy avaliabilities
    for i = 1:length(ηs)
        # And one for each ω parameter
        for j = 1:length(ωs)
            # Can finally generate microbe
            var[i,j] = make_Microbe(MC,γm,Kγ,χl,χu,Pb,d,ϕH,Ωvs[j],fd,ωs[j],R,Reacs,
                       [ηs[i]],[kc],[KS],[kr],n,[1.0],j+(i-1)*length(ηs),PID)
        end
    end
    # Choose sensible initial values
    ϕi = 0.01 # Start at low value
    ai = 1e5
    Ci = 1e3
    Ni = 1e3
    # Choose simulation window
    Tmax = 1e6
    # Preallocate array to store maximum ribosome fractions for each case
    maxϕf = zeros(length(ηs),length(ωs))
    maxϕv = zeros(length(ηs),length(ωs))
    # Loop over all the fixed species
    for i = 1:length(ηs)
        println("Fixed Ω run $i")
        for j = 1:length(ωs)
            # Simulate each population
            C, T = sing_pop(ps,Ni,Ci,ai,ϕi,fix[i,j],Tmax)
            # Extract max ribosome fraction
            maxϕf[i,j] = maximum(C[:,5])
        end
    end
    # Loop over all the variable species
    for i = 1:length(ηs)
        for j = 1:length(ωs)
            # Simulate each population
            C, T = sing_pop(ps,Ni,Ci,ai,ϕi,var[i,j],Tmax)
            # Extract max ribosome fraction
            maxϕv[i,j] = maximum(C[:,5])
        end
    end
    # Check if directory exists and if not make it
    if ~isdir("Output/Fig8")
        mkdir("Output/Fig8")
    end
    # Setup plotting
    pyplot(dpi=200)
    # Load in color scheme
    a = ColorSchemes.sunset.colors
    # Define latex label for max ribosome fraction
    Φm = L"\phi_R"
    ylb = "Maximum $(Φm)"
    # And a latex label for the title
    Ω12 = L"\Omega_{\frac{1}{2}}"
    # EDIT PLOTS TO LOOK NICER
    # First plot fixed case
    p1 = plot(ylim=(0.0,0.5),xlabel="Energy availability",ylabel=ylb)
    plot!(p1,title="Fixed $(Ω12)")
    plot!(p1,ηs,maxϕf[:,1],label="ω = $(ωs[1])",color=a[1])
    plot!(p1,ηs,maxϕf[:,2],label="ω = $(ωs[2])",color=a[2])
    plot!(p1,ηs,maxϕf[:,3],label="ω = $(ωs[3])",color=a[3])
    plot!(p1,ηs,maxϕf[:,4],label="ω = $(ωs[4])",color=a[4])
    plot!(p1,ηs,maxϕf[:,5],label="ω = $(ωs[5])",color=a[5])
    savefig(p1,"Output/Fig8/fixcomp.png")
    # Then plot corresponding variable case
    p2 = plot(ylim=(0.0,0.5),legend=false,xlabel="Energy availability")
    plot!(p2,ylabel=ylb,title="Variable $(Ω12)")
    plot!(p2,ηs,maxϕv[:,1],label="ω = $(ωs[1])",color=a[1])
    plot!(p2,ηs,maxϕv[:,2],label="ω = $(ωs[2])",color=a[2])
    plot!(p2,ηs,maxϕv[:,3],label="ω = $(ωs[3])",color=a[3])
    plot!(p2,ηs,maxϕv[:,4],label="ω = $(ωs[4])",color=a[4])
    plot!(p2,ηs,maxϕv[:,5],label="ω = $(ωs[5])",color=a[5])
    savefig(p2,"Output/Fig8/varcomp.png")
    # Add annotations
    px, py = annpos([0.5;20],[0.0;0.5],0.075,0.0)
    annotate!(p1,px,py,text("A",17,:black))
    px, py = annpos([0.5;20],[0.0;0.5],0.075,0.0)
    annotate!(p2,px,py,text("B",17,:black))
    # Plot all graphs as a single figure
    pt = plot(p1,p2,layout=(2,1),size=(600,800),margin=5.0mm)
    savefig(pt,"Output/Fig8/figure8.png")
    return(nothing)
end

@time figure8()
