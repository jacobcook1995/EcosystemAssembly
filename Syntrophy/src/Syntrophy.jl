module Syntrophy

# Can import or use other modules
# Can also include code from other files other files

# Export objects and function that I want to be externally accessible
export GFree

# Decleration of internally used constants
Rgas = 8.31446261815324 # gas constant in J.K^-1.mol^-1

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
