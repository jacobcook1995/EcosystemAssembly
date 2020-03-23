# This script is similar to MarsParameters but is for the inhibitory model

# export very simple functions
export ↦

# Overloadable version of getfields
↦(s, f) = getfield(s, f)

export InhibParameters, make_InhibParameters, make_Reaction, Reaction, make_Microbe, Microbe

"""
    Reaction(ID::Int64,Rct::Int64,Prd::Int64,ΔG0::Float64)
Type containing the parameters for a particular reaction.
# Arguments
- `ID::Int64`: Number to identify reaction
- `Rct::Int64`: Identity number of reactant
- `Prd::Int64`: Identity number of product
- `ΔG0::Float64`: Standard Gibbs free energy change of the reaction
"""
struct Reaction
    ID::Int64;
    Rct::Int64;
    Prd::Int64;
    ΔG0::Float64;
end

"""
    make_Reaction(ID::Int64,Rct::Int64,Prd::Int64,ΔG0::Float64)
Helper function used internally. Takes values for parameters and returns a `Reaction`object.
Also does checks internally to make sure the values are correct.
"""
function make_Reaction(ID::Int64,Rct::Int64,Prd::Int64,ΔG0::Float64)
    @assert ID > 0 "Reaction must be given a positive ID"
    @assert Rct > 0 "All reactants must have postive IDs"
    @assert Prd > 0 "All products must have postive IDs"
    @assert Prd != Rct "Reactions cannot have same reactant and product"

    return(Reaction(ID,Rct,Prd,ΔG0))
end

"""
    Microbe(m::Float64,g::Float64,R::Int64,Reacs::Vector{Int64},η::Vector{Float64})
Type containing the parameters for a particular microbial strain.
# Arguments
- `m::Float64`: Maintenance energy cost of microbe
- `g::Float64`: Proportionality constant between energy and biomass
- `R::Int64`: Number of reactions
- `Reacs::Vector{Int64}`: Reaction numbers
- `η::Vector{Float64}`: ATP generated per mole of reaction
- `qm::Vector{Float64}`: Maximal rates of substrate consumption
- `KS::Vector{Float64}`: Substrate saturation constants
- `kr::Vector{Float64}`: reversibility factors
"""
struct Microbe
    m::Float64;
    g::Float64;
    R::Int64;
    Reacs::Vector{Int64}
    η::Vector{Float64}
    qm::Vector{Float64}
    KS::Vector{Float64}
    kr::Vector{Float64}
end

"""
    make_Microbe(m::Float64,g::Float64,R::Int64,Reacs::Vector{Int64},η::Vector{Float64},
    qm::Vector{Float64},KS::Vector{Float64},kr::Vector{Float64})
Helper function used internally. Takes values for parameters and returns a `Reaction`object.
Also does checks internally to make sure the values are correct.
"""
function make_Microbe(m::Float64,g::Float64,R::Int64,Reacs::Vector{Int64},η::Vector{Float64},
                        qm::Vector{Float64},KS::Vector{Float64},kr::Vector{Float64})
    # Check that physical parameters have been provided
    @assert m >= 0.0 "Maintenance energy cost (m) cannot be negative"
    @assert g >= 0.0 "Proportionality between energy and biomass (g) cannot be negative"
    @assert R > 0 "Number of reactions must be postive"
    # Check that the vectors have the right length
    @assert length(Reacs) == R "Vector of reactions is the wrong length"
    @assert length(η) == R "Vector of ATP generation rates (η) is the wrong length"
    @assert length(qm) == R "Vector of maximal rates (qm) is the wrong length"
    @assert length(KS) == R "Vector of saturation constants (KS) is the wrong length"
    @assert length(kr) == R "Vector of reversibilities (kr) is the wrong length"
    # Check that values in these vectors are plausible
    @assert all(η .>= 0.0) "η values cannot reasonably be negative"
    @assert all(qm .>= 0.0) "Maximal rates cannot be negative"
    @assert all(KS .>= 0.0) "Saturation constants cannot be negative"
    @assert all(kr .>= 0.0) "Cannot have negative values for the reversibility"
    @assert all(Reacs .> 0) "All reactions are supposed to be indicated by a positive numbers"
    return(Microbe(m,g,R,Reacs,η,qm,KS,kr))
end

"""
    InhibParameters(N::Int64,M::Int64,O::Int64,T::Float64,κ::Vector{Float64},δ::Vector{Float64},
    reacs::Vector{Reaction},mics::Vector{Microbe})
Type containing the parameters for a simuation.
# Arguments
- `N::Int64`: Number of Consumers
- `M::Int64`: Number of Resources
- `O::Int64`: Number of Reactions
- `T::Float64`: Temperature of system
- `κ::Vector{Float64}`: Vector of external resource supply rates
- `δ::Vector{Float64}`: Vector of resource decay rates
- `reacs::Vector{Reaction}`: Vector of reactions
- `mics::Vector{Microbe}`: Vector of microbial strains
"""
struct InhibParameters
    N::Int64;
    M::Int64;
    O::Int64;
    T::Float64;
    κ::Vector{Float64}
    δ::Vector{Float64}
    reacs::Vector{Reaction}
    mics::Vector{Microbe}
end

"""
    make_InhibParameters(N::Int64,M::Int64,O::Int64,T::Float64,κ::Vector{Float64},δ::Vector{Float64},
    reacs::Vector{Reaction},mics::Vector{Microbe})
Helper function used internally. Takes values for parameters and returns a `InhibParameters`object.
Also does checks internally to make sure the values are correct.
"""
function make_InhibParameters(N::Int64,M::Int64,O::Int64,T::Float64,κ::Vector{Float64},δ::Vector{Float64},reacs::Vector{Reaction},mics::Vector{Microbe})
    # Use asserts to ensure that arrays are correct sizes
    @assert length(κ) == M "Vector of external resource supplies (κ) is the wrong length"
    @assert length(δ) == M "Vector of decay rates (δ) is the wrong length"
    @assert length(reacs) == O "Vector of reactions is the wrong length"
    @assert length(mics) == N "Vector of microbes is the wrong length"

    # Use asserts to avoid unphysical values
    @assert all(δ .>= 0) "One or more decay rates in δ are negative"
    @assert T > 0.0 "Temperature must be positive"

    # Check that reaction numbering is reasonable
    @assert ((reacs.↦:ID) == collect(1:O)) "Reaction IDs must be sequential from 1 upwards"
    @assert all((reacs.↦:Rct) .<= M) "Reactant numbers cannot exceed number of metabolites"
    @assert all((reacs.↦:Prd) .<= M) "Product numbers cannot exceed number of metabolites"
    # Check that no reaction is repeated
    for i = 1:O
        @assert all(((reacs.↦:Rct)[i] .!= (reacs.↦:Rct)[1:end .!= i]) .| ((reacs.↦:Prd)[i] .!= (reacs.↦:Prd)[1:end .!= i])) "Reaction $i is repeated"
    end

    # Check that the reactions given in the vector of microbes exist in the vector of reactions
    for i = 1:N
        @assert all((mics.↦:Reacs)[i] .<= O) "Microbe $i assigned to reaction that doesn't exist"
    end
    return(InhibParameters(N,M,O,T,κ,δ,reacs,mics))
end
