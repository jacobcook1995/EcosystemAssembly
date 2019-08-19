using Syntrophy
using Plots
using DifferentialEquations
using LaTeXStrings
import PyPlot

# This is a script to store the functions that plot the figures shown in my Syntrophy SI

# A function to make a plot of the maxcosump rate for a microbial population in glucose
function maxcosump()
    # Nutrient variables
    α = 5.55*10^(-6)
    δ = 2.00*10^(-4) #1.00*10^(-4)
    # make nutrients
    # 1 = glucose, 2 = oxegen, 3 = bicarbonate, 4 = hydrogen ion
    nuts = [Nut(1,false,α,δ),Nut(2,true,0,0),Nut(3,false,0,δ),Nut(4,true,0,0)]
    # Now make reactions
    ΔG0 = -2843800.0
    reac = [React(1,[1,2,3,4],[-1,-6,6,6],ΔG0)]
    # make set of microbes
    ηt = [collect(0:0.5:3.5);collect(4:4:32);collect(36:39);collect(39.25:0.25:43)]
    N = length(ηt) # (maximal η)-1
    r = 1 # Only reaction
    m = 2.16*10^(-19) # maintainance
    mics = Array{Microbe,1}(undef,N)
    for i = 1:N
        mics[i] = Microbe(ηt[i],m,r,0.0)
    end
    # Set intial population and nutrient concentrations for all cases
    pops = 100.0
    concs = zeros(length(nuts))
    # define initial concentrations
    concs[1] = 0.0555 # high initial concentration to ensure growth
    concs[2] = 0.21 # High value so oxegen isn't limiting
    concs[3] = 0.0 # No initial concentration
    concs[4] = 1.00*10.0^(-7) # pH 7
    # Define some constants
    ΔGATP = 75000.0 # Gibbs free energy of formation of ATP in a standard cell
    Temp = 312.0 # Temperature that growth is occuring at in Kelvin
    # All other constants require quite a bit of defining
    Y = 2.36*10.0^(13) # yield in cells per mole of ATP
    KS = 2.40*10.0^(-5) # saturation constant (substrate)
    qm = 3.42*10.0^(-18) # maximal rate substrate consumption mol cell s^-1
    p = [Y,KS,qm,ΔGATP,Temp]
    u0 = [concs;pops] # u0 and p same for every microbe
    tspan = (0.0,5000000.0)
    stoc = (reac.↦:stc)[1]

    # Preallocate vector and calculate maxium ATP generation rate
    atpgen = zeros(N)
    ηs = zeros(N)
    θs = zeros(N)
    for i = 1:N
        settle = false
        # New version of function each time for each microbe
        f(du,u,p,t) = singlepop(du,u,p,nuts,reac,mics[i],t)
        # This step will take a while
        prob = ODEProblem(f,u0,tspan,p)
        sol = solve(prob,adaptive=false,dt=100) # turned dt down to make plots look nicer
        # Check if pop has changed more than 0.1% in last 10 time steps
        diff = (sol'[end,5]-sol'[end-10,5])/sol'[end,5]
        if diff < 0.001
            settle = true
        end
        # Print if it hasn't settled down to sufficently steady popultaion
        if settle == false
            println("Problem with eta = $(ηt[i])")
            println(diff)
        end
        # Now calulate max ATP genration rate
        ηc = mics[i].η # ATP usage amount
        qr = qrate(sol'[end,1:4],KS,qm,ΔGATP,ΔG0,Temp,stoc,ηc) # reaction rate
        atpgen[i] = sol'[end,5]*ηc*qr # multiply by population to get total rate
        θs[i] = θT(sol'[end,1:4],stoc,ΔGATP,ΔG0,ηc,Temp) # obtain thermodynamic term
        ηs[i] = ηc # save η for later plotting
    end
    pyplot(dpi=150)
    # Plot maximal ATP generation rate
    plot(ηs,atpgen,xaxis=L"\eta\;\;mol_{ATP}\;(mol_{reaction})^{-1}",yaxis="Maximal ATP production rate mol/s",label="")
    plot!(title="ATP rate vs efficency trade off")
    savefig("Output/Atpgenrate.png")
    # Plot thermodynamicinhibition term
    plot(ηs,θs,xaxis=L"\eta\;\;mol_{ATP}\;(mol_{reaction})^{-1}",yaxis=L"\theta\;\;",label="")
    plot!(title="Thermodynamic inhibition as a function of efficency")
    savefig("Output/ThermInhib.png")
    return(nothing)
end

@time maxcosump()
