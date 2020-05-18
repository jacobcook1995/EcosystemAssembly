# Script that makes the parameters needed for simulation of the single species proteome model
# Probably should get merged eventually

export initialise_prot

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
    # Now make the parameter set
    ps = make_ProtParameters(MC,γm,T,n)
    return(ps)
end
