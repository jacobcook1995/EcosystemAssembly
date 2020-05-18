# Script that makes the parameters needed for simulation of the single species proteome model
# Probably should get merged eventually

export initialise_prot, make_var_prot

# function to generate parameter set for the model with inhibition
function initialise_prot()
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # Cell mass is taken from Bremer H, Dennis P (1996) Modulation of chemical
    # composition and other parameters of the cell by growth rate (Book chapter).
    MC = 10^8
    # Elongation rate also taken from Bremer (1996)
    γm = 1260.0
    # Now preallocate protein masses
    n = zeros(Int64,3)
    # ribosome mass taken from Keseler IM, et al. (2011)
    n[1] = 7459
    # Other protein mass averaged from Brandt F, et al. (2009)
    n[2:3] .= 300
    # Now make the reaction
    ΔG = -6e5 # Relatively small Gibbs free energy change
    r = make_Reaction(1,1,2,ΔG)
    # η chosen so that a substantial portion of the Gibbs free energy is retained
    η = 0.9*(-ΔG/ΔGATP)
    # Now make the parameter set
    ps = make_ProtParameters(MC,γm,T,η,r,n)
    return(ps)
end

# A function to make set of parameters that are varied in our simulation
function make_var_prot(ps::ProtParameters,ϕ::Array{Float64,1})
    # Check that housekeeping fraction matches with Scott et al. 2010
    @assert ϕ[3] == 0.45 "Housekeeping protein fraction should be 0.45"
    # Find number of enzymes using the protein fraction
    E = Eα(ϕ[2],ps)
    # Now make parameter set
    pa = make_VarProtParameters(ϕ,E)
    return(pa)
end
