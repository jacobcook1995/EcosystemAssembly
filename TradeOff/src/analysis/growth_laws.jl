using TradeOff
using Plots
using LaTeXStrings
using Random
using LsqFit
import PyPlot

function new_mic_gl(γm::Float64,ω::Float64,χ::Float64,Kγ::Float64,KΩ::Float64,M::Int64,ps::TOParameters)
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
    # Always want to allow near to equilbrium reactions now
    syn = true
    # Use formula to calculate how many reactions are implied
    O = 2*M - 3
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # Cell mass is taken from Bremer H, Dennis P (1996) Modulation of chemical
    # composition and other parameters of the cell by growth rate (Book chapter).
    MC = 10^8
    # Now preallocate protein masses
    n = zeros(Int64,3)
    # ribosome mass taken from Keseler IM, et al. (2011)
    n[1] = 7459
    # Other protein mass averaged from Brandt F, et al. (2009)
    n[2:3] .= 300
    # Choosing a low value to better highlight the growth laws
    d = 6.0e-10
    # The proportion of ribosomes bound is taken from Underwood et al to be 70%
    Pb = 0.7
    # Housekeeping fraction is taken from Scott et al. 2010
    ϕH = 0.45
    # Number of doublings required to dilute to 1%
    fd = log(100)/log(2)
    # Each strain should just have the one available reaction
    R = 1
    Reacs = [1]
    # Make vectors of the (fixed) kinetic parameters
    kcs = kc.*ones(R)
    KSs = KS.*ones(R)
    krs = kr.*ones(R)
    # Reactions given random proportional weightings, done this in the simplest way possible
    ϕP = rand(R)
    ϕP = ϕP/sum(ϕP)
    # Find corresponding η's for these reactions
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
    mic = new_make_Microbe(MC,γm,Kγ,χ,Pb,d,ϕH,KΩ,fd,ω,R,Reacs,η,kcs,KSs,krs,n,ϕP,ID,PID)
    return(mic)
end

# function to generate parameter set for the fixed parameters
function initialise_gl(M::Int64,O::Int64,μrange::Float64)
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # Want waste to be removed but not the substrate
    δ = [0.0,6.0e-5]
    # Human blood glucose is approximatly 5.5 m mol per litre (wikipedia)
    # Sensible order of magnitude to aim for set κ/δ = 5.5e-3
    κ = zeros(M)
    # All but resource 1 is not supplied
    κ[1] = 3.3e-7 # Metabolite supply rate
    # Generate fixed set of reactions
    RP, ΔG = fix_reactions(O,M,μrange,T)
    # Preallocate vector of reactions
    reacs = Array{Reaction,1}(undef,O)
    for i = 1:O
        reacs[i] = make_Reaction(i,RP[i,1],RP[i,2],ΔG[i])
    end
    # Now make the parameter set
    ps = make_TOParameters(M,O,T,κ,δ,reacs)
    return(ps)
end

# function to plot growth laws
function growth_laws()
    println("Successfully compiled.")
    # Simple test data set
    ai = 5.0 # initial energy level
    Ni = 100.0 # initial population
    Si = 1.0
    Pi = 0.0
    ϕi = 0.1
    # These are parameters I might need to change to get reasonable results
    ω = 1.0
    χ = 30.0
    # These require real thought
    Kγ = 1.0e12
    KΩ = 5.0e10
    # Simple 1 reaction case
    M = 2
    O = 1
    # Choose simulation time
    Tmax = 5e6
    # Set this as a middling value of Gibbs free energy change
    μr = 3e6 # Relatively small Gibbs free energy change
    # Max Elongation rate also taken from Bremer (1996), convert from minutes to seconds
    γm = 1260.0/60.0 # Change this one
    # Make vector of γm values
    γms = [γm/100,γm/50,γm/20,γm/10,γm/7.5,γm/5,γm/4,γm/3,γm/2,γm/1.5,γm]
    # Make vector to store final growth rates and fractions
    λ1 = zeros(length(γms))
    ϕ1 = zeros(length(γms))
    # Setup plotting options
    pyplot(dpi=200)
    pR = L"\phi_R"
    p1 = plot(xaxis="Growth rate, λ",yaxis="Ribosome fraction, $(pR)")
    for i = 1:length(γms)
        # Initialise parameters
        ps = initialise_gl(M,O,μr)
        # and then the microbe
        mic = new_mic_gl(γms[i],ω,χ,Kγ,KΩ,M,ps)
        # Simulate dynamics of the single population
        C, T = sing_pop(ps,Ni,Si,ai,ϕi,mic,Tmax)
        # Now calculate growth rates and elongation rates
        λa = zeros(length(T))
        for j = 1:length(T)
            λa[j] = λs(C[j,4],C[j,5],mic)
        end
        # Save maximum λ and ϕR values
        λ1[i], ind = findmax(λa)
        ϕ1[i] = C[ind,5]
    end
    # Now make set of ΔG values
    μrs = [μr,μr/1.5,μr/2,μr/3,μr/4,μr/5,μr/7.5,μr/10,μr/20,μr/50,μr/100]
    # Make vector to store final growth rates and fractions
    λ2 = zeros(length(μrs))
    ϕ2 = zeros(length(μrs))
    for i = 1:length(μrs)
        # Initialise parameters
        ps = initialise_gl(M,O,μrs[i])
        # and then the microbe
        mic = new_mic_gl(γm,ω,χ,Kγ,KΩ,M,ps)
        # Simulate dynamics of the single population
        C, T = sing_pop(ps,Ni,Si,ai,ϕi,mic,Tmax)
        # Now calculate growth rates and proteome fractions
        λa = zeros(length(T))
        for j = 1:length(T)
            λa[j] = λs(C[j,4],C[j,5],mic)
        end
        # Save final λ and ϕR values
        λ2[i], ind = findmax(λa)
        ϕ2[i] = C[ind,5]
    end
    # Now want to do a least squares fit for both sets of data
    @. model(x, p) = p[1]*x + p[2]
    p0 = [0.5,0.5]
    fit1 = curve_fit(model,λ1[4:end],ϕ1[4:end],p0)
    pr1 = coef(fit1)
    fit2 = curve_fit(model,λ2,ϕ2,p0)
    pr2 = coef(fit2)
    # plot both lines on the graph
    λ1s = [0.0;λ1]
    λ2s = [0.0;λ2]
    plot!(p1,λ1s,pr1[1]*λ1s.+pr1[2],label="",color=:blue)
    plot!(p1,λ2s,pr2[1]*λ2s.+pr2[2],label="",color=:red)
    # Add arrows indicating direction of change
    l = 1e-4 # Way too large
    quiver!(p1,[λ1[end-2]],[ϕ1[end-2]+0.03],quiver=([-l],[-pr1[1]*l]),color=:blue)
    quiver!(p1,[λ2[6]],[ϕ2[6]-0.03],quiver=([l],[pr2[1]*l]),color=:red)
    # Position is where the annotation centres are
    pos1x = λ1[end-2] - l/2
    pos1y = ϕ1[end-2] + 0.05 - pr1[1]*l/2
    pos2x = λ2[6] + l/2
    pos2y = ϕ2[6] - 0.05 + pr2[1]*l/2
    # Calculate rotations in degrees
    r1 = -12.5 # Needs to be -ve
    r2 = 20 # Needs to be +ve
    # Then add the annotations
    annotate!(p1,pos1x,pos1y,text("Translational inhibition",8,color=:blue,rotation=r1))
    annotate!(p1,pos2x,pos2y,text("Nutrient quality",8,color=:red,rotation=r2))
    # Plot final values
    scatter!(p1,λ1,ϕ1,label="",color=:lightblue,markersize=5)
    scatter!(p1,λ2,ϕ2,label="",color=:orange,markersize=5)
    # Finally save the graph
    savefig(p1,"Output/GrowthLaws.png")
    return(nothing)
end

@time growth_laws()
