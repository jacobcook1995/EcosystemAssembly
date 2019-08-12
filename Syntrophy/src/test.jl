using Syntrophy
using Plots

# This is a script to write my testing code into
# Anything reusable should be moved into Syntrophy.jl as a seperate callable function

# function to calculate the rate of substrate consumption q
function qrate(concs::Array{Float64,1},K::Float64,qm::Float64,ΔGATP::Float64,
                    ΔG0::Float64,Temp::Float64,stoc::Array{Int64,1},η::Float64)
    # concs => Vector of nutrient concentrations
    # K => Saturation constant for the substrate
    # qm => Maximal reaction rate for substrate
    # stoc => stochiometry vector
    # ΔGATP => Gibbs free energy to form ATP in standard cell
    # ΔG0 => Standard gibbs free energy of reaction
    # Temp => Temperature in Kelvin
    # η => free energy use strategy, mol of ATP per mol of substrate
    ############ START OF FUNCTION ###################

    # calulate substrate coefficent
    S = SubCoef(concs,stoc)
    # Call function to find thermodynamic factor θ
    θ = θT(concs,stoc,ΔGATP,ΔG0,η,Temp)
    # Only η changes between species
    q = qm*S*(1-θ)/(K+S*(1+θ))
    return(q)
end



function gluc()
    # Nutrient variables
    α = 5.55*10^(-6)
    δ = 1.00*10^(-4)
    # make nutrients
    # 1 = glucose, 2 = oxegen, 3 = bicarbonate, 4 = hydrogen ion
    nuts = [Nut(1,false,α,δ),Nut(2,true,0,0),Nut(3,false,0,δ),Nut(4,true,0,0)]
    # Now make reactions
    ΔG0 = -2843800.0
    reac = [React(1,[1,2,3,4],[-1,-6,6,6],ΔG0)]
    # microbe variables
    η = 32.0
    r = 1 # Only reaction
    # Considering 1 microbe with no maintaince and no dilution
    mics = [Microbe(η,0.0,1,0.0)]
    # Set intial populations and nutrient concentrations
    pops = 10.0*ones(length(mics))
    concs = zeros(length(nuts))
    # define initial concentrations
    concs[1] = 0.0555 # high initial concentration to ensure growth
    concs[2] = 0.0375 # High value so oxegen isn't limiting
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
    u0 = [concs;pops]

    # Testing qrate
    stoc = (reac.↦:stc)[1]
    q = qrate(concs,KS,qm,ΔGATP,ΔG0,Temp,stoc,η)
    return(nothing)
end

@time gluc()
