using Assembly
using Plots
using Plots.PlotMeasures
import PyPlot

# Function to visualise syntrophy and pollution interactions
function int_vis()
    println("Successfully compiled")
    # Assume that half saturation occurs at a quarter κ/δ
    KS = (1/4)*5.5e-3
    # Arbitary number that seems to give decent survival
    kc = 10.0
    # The reversibility factor remains the same as previously
    kr = 10.0
    # Initial ribosome fraction is taken from my ATP fits
    ϕR0 = 0.128
    # Initial energy concentration
    a0 = 1e5
    # Make relatively short Tmax
    Tmax = 5e5
    # Make parameter set for syntrophy
    ps = initialise_syn_trio(kc,KS,kr)
    # Initial conditions for no sytrophy case
    pop = [1e3,1e3,0.0]
    conc = [2.5e-3,2.5e-3,0.0]
    as = a0*ones(3)
    ϕs = ϕR0*ones(3)
    # Then run this simulation
    Cns, Tns = full_simulate(ps,Tmax,pop,conc,as,ϕs)
    # Initial conditions for with sytrophy case
    pop = [1e3,1e3,1e3]
    # Then run this simulation
    Cws, Tws = full_simulate(ps,Tmax,pop,conc,as,ϕs)
    # Make parameter set for syntrophy
    ps = initialise_pol_duo(kc,KS,kr)
    # Initial conditions no pollution case
    pop = [1e3,0.0]
    as = a0*ones(2)
    ϕs = ϕR0*ones(2)
    # Then run this simulation
    Cnp, Tnp = full_simulate(ps,Tmax,pop,conc,as,ϕs)
    # Initial conditions with pollution case
    pop = [1e3,1e3]
    # Then run this simulation
    Cwp, Twp = full_simulate(ps,Tmax,pop,conc,as,ϕs)
    # initialise the plotting
    pyplot(dpi=200)
    # First population plot
    p1 = plot(xlabel="Time (s)",ylabel="Population",yaxis=:log10,legend=:bottomright)
    plot!(p1,Tns,Cns[:,1:2],labels=["Species A" "Species B"],title="No syntrophy interaction")
    # First concentration plot
    p2 = plot(xlabel="Time (s)",ylabel="Concentration",legend=:right)
    plot!(p2,Tns,Cns[:,4:6],labels=["Metabolite 1" "Metabolite 2" "Metabolite 3"])
    # Second population plot
    p3 = plot(xlabel="Time (s)",ylabel="Population",yaxis=:log10,legend=:bottomright)
    plot!(p3,Tws,Cws[:,1:3],labels=["Species A" "Species B" "Species C"],title="Syntrophy interaction")
    # Second concentration plot
    p4 = plot(xlabel="Time (s)",ylabel="Concentration",legend=:right)
    plot!(p4,Tws,Cws[:,4:6],labels=["Metabolite 1" "Metabolite 2" "Metabolite 3"])
    # Combine all the plots
    pc = plot(p1,p2,p3,p4,layout=4,size=(1200,800),margin=7.5mm)
    savefig(pc,"Output/SI/SynPlot.png")
    # First population plot
    p1 = plot(xlabel="Time (s)",ylabel="Population",yaxis=:log10,legend=:bottomright)
    plot!(p1,Tnp,Cnp[:,1],label="Species A",title="No pollution interaction")
    # First concentration plot
    p2 = plot(xlabel="Time (s)",ylabel="Concentration",legend=:right)
    plot!(p2,Tnp,Cnp[:,3:5],labels=["Metabolite 1" "Metabolite 2" "Metabolite 3"])
    # Second population plot
    p3 = plot(xlabel="Time (s)",ylabel="Population",yaxis=:log10,legend=:bottomright)
    plot!(p3,Twp,Cwp[:,1:2],labels=["Species A" "Species B"],title="Pollution interaction")
    # Second concentration plot
    p4 = plot(xlabel="Time (s)",ylabel="Concentration",legend=:right)
    plot!(p4,Twp,Cwp[:,3:5],labels=["Metabolite 1" "Metabolite 2" "Metabolite 3"])
    # Combine all the plots
    pc = plot(p1,p2,p3,p4,layout=4,size=(1200,800),margin=7.5mm)
    savefig(pc,"Output/SI/PolPlot.png")
    return(nothing)
end

@time int_vis()
