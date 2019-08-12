module Syntrophy

# Can import or use other modules
# Can also include code from other files other files

# Export complicated functions
export GFree, θT, netE
# expoert all my objects
export Nut, React, Microbe
# export very simple functions
export ↦

# Decleration of internally used constants
Rgas = 8.31446261815324 # gas constant in J.K^-1.mol^-1

# Structure to store nutrient identies, supply rates and washout fractions
struct Nut
    idt::Int64 # Identifying integer should be unique
    cst::Bool # Bool to set concentration constant or not
    α::Float64 # Supply rate mol l^−1 s^−1
    δ::Float64 # dilution rate s^-1
    # Test to ensure that constant nutrients don't have supply or dilution rates
    Nut(x,y,z,h) = (y == true && (z != 0.0 || h != 0.0)) ? error("Constant nutrient?") : new(x,y,z,h)
end

# Structure to store details of each reaction
struct React
    idt::Int64 # Identifying integer should be unique
    nidt::Array{Int64,1} # nutrient identifiers
    stc::Array{Int64,1} # stochiometry -ve indicates substrate
    ΔG0::Float64 # Gibbs free energy (in Joules) at standard conditions
end

# Structure to store microbial strategies and populations
struct Microbe
    η::Float64 # free energy use strategy, mol of ATP per mol of substrate
    m::Float64 # maintainance cost, mol of ATP per cell per second
    reac::Int64 # Identity of reaction used to fuel
    δ::Float64 # dilution rate s^-1
end

# Overloadable version of getfields
↦(s, f) = getfield(s, f)

# Function to calculate Gibbs free energy change (in Joules) of the reaction
function GFree(concs::Array{Float64,1},stoc::Array{Int64,1},Temp::Float64,ΔG0::Float64)
    # concs => concentrations in Moles
    # stoc => reaction stochiometries
    # Temp => temperature
    # ΔG0 => Gibbs free energy change in standard conditions
    ############ START OF FUNCTION ###################

    # Expecting stochiometry of one reaction and corresponding concentrations
    # so should be vectors of same length
    if length(stoc) != length(concs)
        error("Data of mismatching length!")
    end

    # Find reaction quotient, Q
    Q = 1
    for i = length(stoc)
        if stoc[i] > 0
            Q *= concs[i]^(stoc[i])
        elseif stoc[i] < 0
            Q *= concs[i]^(stoc[i])
        else
            error("Should not provide any non reacting species in reaction")
        end
    end

    # Calculate temp dependant factor
    RT = Rgas*Temp

    # Find and return Gibbs free energy change
    ΔG = ΔG0 + RT*log(Q)
    return(ΔG)
end

# function to calulate thermodynamic term θ
function θT(concs::Array{Float64,1},stoc::Array{Int64,1},ΔGATP::Float64,ΔG0::Float64,η::Float64,Temp::Float64)
    # concs => Vector of nutrient concentrations
    # stoc => stociometry vector
    # ΔGATP => Gibbs free energy of ATP in a cell
    # ΔG0 => standard Gibbs free energy of the reaction
    # η => free energy use strategy, mol of ATP per mol of substrate
    # Temp => Temperature in Kelvin
    ############ START OF FUNCTION ###################

    # Expecting stochiometry of one reaction and corresponding concentrations
    # so should be vectors of same length
    if length(stoc) != length(concs)
        error("Data of mismatching length!")
    end

    # Find reaction quotient, Q
    Q = 1
    for i = length(stoc)
        if stoc[i] > 0
            Q *= concs[i]^(stoc[i])
        elseif stoc[i] < 0
            Q *= concs[i]^(stoc[i])
        else
            error("Should not provide any non reacting species in reaction")
        end
    end

    # Calculate temp dependant factor
    RT = Rgas*Temp

    θ = Q*exp((ΔG0+η*ΔGATP)/RT)
    return(θ)
end

# Function to calulatenet energy retained by cell after maintainance contribution
function netE(η::Float64,rate::Float64,m::Float64)
    # rate => rate of reaction
    # m => # maintainance cost, mol of ATP per cell per second
    # η => free energy use strategy, mol of ATP per mol of substrate
    ############ START OF FUNCTION ###################

    E = η*rate - m
    return(E)
end

end # module
