module Syntrophy

# Can import or use other modules
# Can also include code from other files other files
using DifferentialEquations

# Export objects and function that I want to be externally accessible
export GFree, Nut, React, Microbe, ↦

# Decleration of internally used constants
Rgas = 8.31446261815324 # gas constant in J.K^-1.mol^-1

# Structure to store nutrient identies, supply rates and washout fractions
struct Nut
    idt::Int64 # Identifying integer should be unique
    α::Float64 # Supply rate mol ml^−1 s^−1
    δ::Float64 # dilution rate s^-1
end

# Structure to store details of each reaction
struct React
    idt::Int64 # Identifying integer should be unique
    sub::Array{Int64,1} # substrate identifier
    subS::Array{Int64,1} # substrate stochiometry
    prd::Array{Int64,1} # end product identifier
    prdS::Array{Int64,1} # product stochiometry
    ΔG0::Float64 # Gibbs free energy (in Joules) at standard condidtions
end

# Structure to store microbial strategies and populations
struct Microbe
    η::Float64 # free energy use strategy, mol of ATP per mol of substrate
    m::Float64 # maintainance cost, mol of ATP per cell per second
    reac::Int64 # Identity of reaction used to fuel
end

# Overloadable version of getfields
↦(s, f) = getfield(s, f)

# Function to calculate Gibbs free energy change (in Joules) of the reaction
function GFree(subC::Array{Float64,1},subN::Array{Int64,1},prodC::Array{Float64,1},prodN::Array{Int64,1},Temp::Float64,ΔG0::Float64)
    # SubC => substrate concentrations in Moles
    # SubN => substrate stochiometry
    # ProdC => product concentrations in Moles
    # ProdN => product stochiometry
    # Temp => temperature
    # ΔG0 => Gibbs free energy change in standard conditions
    ############ START OF FUNCTION ###################

    # Calculate reaction quotient, Q
    den = 1.0
    for i = 1:length(subC)
        den *= (subC[i])^subN[i]
    end
    num = 1.0
    for i = 1:length(prodC)
        num *= (prodC[i])^prodN[i]
    end
    Q = num/den

    # Calculate temp dependant factor
    RT = Rgas*Temp

    # Find and return Gibbs free energy change
    ΔG = ΔG0 + RT*log(Q)
    return(ΔG)
end

end # module
