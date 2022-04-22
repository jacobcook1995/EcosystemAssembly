# This script of the parameters for the full model, including proteome parameters

# export very simple functions
 export ↦

 # Overloadable version of getfields
 ↦(s, f) = getfield(s, f)

export Reaction, make_Reaction, Microbe, make_Microbe

export TOParameters, make_TOParameters, MicData, make_MicData

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
    Microbe(MC::Int64,γm::Float64,Kγ::Float64,χl::Float64,χu::Float64,Pb::Float64,d::Float64,ϕH::Float64,
    KΩ::Float64,fd::Float64,R::Int64,Reacs::Vector{Int64},η::Vector{Float64},kc::Vector{Float64},
    KS::Vector{Float64},kr::Vector{Float64},n::Vector{Int64},ϕP::Vector{Float64},ID::Int64,PID::String)
Type containing the parameters for a particular microbial strain, this now includes proteome parameters.
# Arguments
- `MC::Int64`: Mass of cell in amino acids.
- `γm::Float64`: Maximum elongation rate in amino acids per minute (per ribosome)
- `Kγ::Float64`: Threshold energy in ATP per cell
- `χl::Float64`: minimum ATP use per translation step
- `χu::Float64`: maximum (additional) ATP use per translation step
- `Pb::Float64`: Proportion of ribosomes bound
- `d::Float64`: Biomass loss rate (~death rate)
- `ϕH::Float64`: Fraction of proteome allocated to housekeeping proteins
- `KΩ::Float64`: Saturation constant for proteome fraction with energy
- `fd::Float64`: Number of doublings required to reorder proteome
- `ω::Float64`: Fraction of maximum ribsome fraction can reach
- `R::Int64`: Number of reactions
- `Reacs::Vector{Int64}`: Reaction numbers
- `η::Vector{Float64}`: ATP generated per mole of reaction
- `kc::Vector{Float64}`: Catalytic rate constants
- `KS::Vector{Float64}`: Substrate saturation constants
- `kr::Vector{Float64}`: reversibility factors
- `n::Vector{Int64}`: Array of number of amino acids per protein type, [r,p,h]
- `ϕP::Vector{Float64}`: Vector of metabolic fractions
- `ID::Int64`: Identifer of microbe in the pool
- `PID::String`: Identifier of the pool
"""
struct Microbe
    MC::Int64;
    γm::Float64;
    Kγ::Float64;
    χl::Float64;
    χu::Float64;
    Pb::Float64;
    d::Float64;
    ϕH::Float64;
    KΩ::Float64;
    fd::Float64;
    ω::Float64;
    R::Int64;
    Reacs::Vector{Int64}
    η::Vector{Float64}
    kc::Vector{Float64}
    KS::Vector{Float64}
    kr::Vector{Float64}
    n::Vector{Int64}
    ϕP::Vector{Float64}
    ID::Int64
    PID::String
end

"""
    make_Microbe(MC::Int64,γm::Float64,Kγ::Float64,χl::Float64,χu::Float64,Pb::Float64,d::Float64,ϕH::Float64,
    KΩ::Float64,fd::Float64,ω::Float64,R::Int64,Reacs::Vector{Int64},η::Vector{Float64},kc::Vector{Float64},
    KS::Vector{Float64},kr::Vector{Float64},n::Vector{Int64},ϕP::Vector{Float64},ID::Int64,PID::String)
Helper function used internally. Takes values for parameters and returns a `Microbe` object.
Also does checks internally to make sure the values are correct.
"""
function make_Microbe(MC::Int64,γm::Float64,Kγ::Float64,χl::Float64,χu::Float64,Pb::Float64,d::Float64,ϕH::Float64,
            KΩ::Float64,fd::Float64,ω::Float64,R::Int64,Reacs::Vector{Int64},η::Vector{Float64},kc::Vector{Float64},
            KS::Vector{Float64},kr::Vector{Float64},n::Vector{Int64},ϕP::Vector{Float64},ID::Int64,PID::String)
    # Check that physical parameters have been provided
    @assert R > 0 "number of reactions must be postive"
    @assert MC > 0 "cell mass must be positive"
    @assert γm >= 0.0 "elongation rate cannot be negative"
    @assert χl > 0.0 "minimum ATP per sythesis step must be positive"
    @assert χu > 0.0 "maximum (additional) ATP per sythesis step must be positive"
    @assert Kγ > 0.0 "threshold energy must be positive"
    @assert d > 0.0 "death rate must be positive"
    @assert 0.0 <= Pb <= 1.0 "proportion of ribosomes bound has to be between 0 and 1"
    @assert 0.0 <= ϕH <= 1.0 "housekeeping proteins have to be between 0 and 100% of total"
    @assert KΩ > 0.0 "proteome fraction saturation constant must be positive"
    @assert fd > 0.0 "number of doublings required must be positive"
    @assert 0.0 <= ω <= 1.0 "ribosome fraction maximum has to be between 0.0 and 1.0 of total"
    # Check that the vectors have the right length
    @assert length(Reacs) == R "Vector of reactions is the wrong length"
    @assert length(η) == R "Vector of ATP generation rates (η) is the wrong length"
    @assert length(kc) == R "Vector of maximal rates (qm) is the wrong length"
    @assert length(KS) == R "Vector of saturation constants (KS) is the wrong length"
    @assert length(kr) == R "Vector of reversibilities (kr) is the wrong length"
    @assert length(ϕP) == R "Vector of metabolic fractions is the wrong length"
    # Check that values in these vectors are plausible
    @assert all(η .>= 0.0) "η values cannot reasonably be negative"
    @assert all(kc .>= 0.0) "Catalytic rate constants cannot be negative"
    @assert all(KS .>= 0.0) "Saturation constants cannot be negative"
    @assert all(kr .>= 0.0) "Cannot have negative values for the reversibility"
    @assert all(Reacs .> 0) "All reactions are supposed to be indicated by a positive numbers"
    @assert length(n) == 2+R "Only considering $(2+R) protein types"
    @assert all(n .> 0) "All proteins must have positive mass"
    @assert all(ϕP .>= 0.0) "All metabolic fractions must be non-negative"
    @assert sum(ϕP) ≈ 1.0 "Metabolic fractions should sum to 1"
    return(Microbe(MC,γm,Kγ,χl,χu,Pb,d,ϕH,KΩ,fd,ω,R,Reacs,η,kc,KS,kr,n,ϕP,ID,PID))
end

"""
    TOParameters(M::Int64,O::Int64,T::Float64,κ::Vector{Float64},δ::Vector{Float64},
    reacs::Vector{Reaction})
Type containing the parameters for a simuation.
# Arguments
- `M::Int64`: Number of Resources
- `O::Int64`: Number of Reactions
- `T::Float64`: Temperature of system
- `κ::Vector{Float64}`: Vector of external resource supply rates
- `δ::Vector{Float64}`: Vector of resource decay rates
- `reacs::Vector{Reaction}`: Vector of reactions
"""
struct TOParameters
    M::Int64;
    O::Int64;
    T::Float64;
    κ::Vector{Float64}
    δ::Vector{Float64}
    reacs::Vector{Reaction}
end

"""
    make_TOParameters(M::Int64,O::Int64,T::Float64,κ::Vector{Float64},δ::Vector{Float64},reacs::Vector{Reaction})
Helper function used internally. Takes values for parameters and returns a `TOParameters`object.
Also does checks internally to make sure the values are correct.
"""
function make_TOParameters(M::Int64,O::Int64,T::Float64,κ::Vector{Float64},δ::Vector{Float64},reacs::Vector{Reaction})
    # Use asserts to ensure that arrays are correct sizes
    @assert length(κ) == M "Vector of external resource supplies (κ) is the wrong length"
    @assert length(δ) == M "Vector of decay rates (δ) is the wrong length"
    @assert length(reacs) == O "Vector of reactions is the wrong length"

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

    return(TOParameters(M,O,T,κ,δ,reacs))
end

"""
    MicData(MID::Int64,PID::String,ImT::Float64,ExT::Float64)
Type containing the parameters for a simuation.
# Arguments
- `MID::Int64`: Microbe ID
- `PID::String`: Pool ID
- `ImT::Float64`: Time of immigration
- `ExT::Float64`: Time of extinction
"""
struct MicData
    MID::Int64;
    PID::String;
    ImT::Float64;
    ExT::Float64;
end

"""
    make_MicData(MID::Int64,PID::String,ImT::Float64,ExT::Float64)
Helper function used internally. Takes values for parameters and returns a `MicData`object.
Also does checks internally to make sure the values are correct.
"""
function make_MicData(MID::Int64,PID::String,ImT::Float64,ExT::Float64)
    # Use asserts to ensure that arrays are correct sizes
    @assert MID > 0 "microbe IDs must be positive"
    @assert ImT >= 0.0 "negative immigration times are not allowed"
    @assert isnan(ExT) || ExT >= ImT "strain cannot go extinct before it arrives"
    @assert length(PID) == 8 "pool ID is supposed to be an 8 character string"

    return(MicData(MID,PID,ImT,ExT))
end
