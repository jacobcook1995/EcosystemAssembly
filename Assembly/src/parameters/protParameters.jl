# This script of the parameters for the single population proteome model
# This should be later integrated with inhibParameters.

export ProtParameters, make_ProtParameters, VarProtParameters, make_VarProtParameters

"""
    ProtParameters(MC::Int64,γm::Float64,T::Float64,η::Float64,KS::Float64,kr::Float64,
    kc::Float64,ρ::Float64,Kγ::Float64,d::Float64,r::Reaction,n::Array{Int64,1},
    δ::Array{Float64,1},κ::Array{Float64,1})
Type containing the parameters for a simulation of the single population proteome model.
# Arguments
- `MC::Int64`: Mass of cell in amino acids.
- `γm::Float64`: Maximum elongation rate in amino acids per minute (per ribosome)
- `T::Float64`: Temperature of the ecosystem
- `η::Float64`: ATP generated per mole of reaction
- `KS::Float64`: Substrate saturation constant
- `kr::Float64`: Reversibility factor
- `kc::Float64`: Catalytic rate constant
- `ρ::Float64`: ATP per sythesis step
- `Kγ::Float64`: Threshold energy in ATP per cell
- `d::Float64`: Biomass loss rate (~death rate)
- `r::Reaction`: Reaction that population uses
- `n::Array{Int64,1}`: Array of number of amino acids per protein type, [r,p,h]
- `δ::Array{Float64,1}`: Dilution rates of the metabolites
- `κ::Array{Float64,1}`: Supply rates of the metabolites
"""
struct ProtParameters
    MC::Int64;
    γm::Float64;
    T::Float64;
    η::Float64;
    KS::Float64;
    kr::Float64;
    kc::Float64;
    ρ::Float64;
    Kγ::Float64;
    d::Float64;
    r::Reaction;
    n::Array{Int64,1}
    δ::Array{Float64,1}
    κ::Array{Float64,1}
end

"""
    make_ProtParameters(MC::Int64,γm::Float64,T::Float64,η::Float64,KS::Float64,kr::Float64,
    kc::Float64,ρ::Float64,Kγ::Float64,d::Float64,r::Reaction,n::Array{Int64,1},δ::Array{Float64,1},
    κ::Array{Float64,1})
Helper function used internally. Takes values for parameters and returns a `ProtParameters`object.
Also does checks internally to make sure the values are correct.
"""
function make_ProtParameters(MC::Int64,γm::Float64,T::Float64,η::Float64,KS::Float64,kr::Float64,
                            kc::Float64,ρ::Float64,Kγ::Float64,d::Float64,r::Reaction,n::Array{Int64,1},
                            δ::Array{Float64,1},κ::Array{Float64,1})
    # Use asserts to ensure parameters a sensible
    @assert MC > 0 "Cell mass must be positive"
    @assert γm >= 0.0 "Elongation rate cannot be negative"
    @assert T > 0.0 "Temperature must be positive"
    @assert η > 0.0 "Cannot generate a negative amount of ATP"
    @assert KS > 0.0 "Negative saturation constants arn't possible"
    @assert kr > 0.0 "Cannot have a negative reversibility factor"
    @assert kc > 0.0 "Catalytic rate constant cannot be negative"
    @assert ρ > 0.0 "ATP per sythesis step must be positive"
    @assert Kγ > 0.0 "Threshold energy must be positive"
    @assert d > 0.0 "Death rate must be positive"

    # Similar asserts for the vector parameters
    @assert length(n) == 3 "Only considering 3 protein types"
    @assert all(n .> 0) "All proteins must have positive mass"
    @assert r.Rct == 1 && r.Prd == 2 "Reaction must be between metabolite 1 and 2"
    @assert length(δ) == 2 && length(κ) == 2 "Only two metabolites"
    @assert all(κ .>= 0) "Negative metabolite supply isn't possible"
    @assert all(δ .> 0) "All dilution rates must be positive"

    return(ProtParameters(MC,γm,T,η,KS,kr,kc,ρ,Kγ,d,r,n,δ,κ))
end

"""
    VarProtParameters(ϕ::Array{Float64,1},E::Float64)
Type containing the parameters varied in a simulation of the single population proteome model.
# Arguments
- `ϕ::Array{Float64,1}`: Protein fractions.
- `E::Float64`: Number of enzymes per cell
"""
struct VarProtParameters
    ϕ::Array{Float64,1}
    E::Float64;
end

"""
    make_VarProtParameters(ϕ::Array{Float64,1},E::Float64)
Helper function used internally. Takes values for parameters and returns a `VarProtParameters`object.
Also does checks internally to make sure the values are correct.
"""
function make_VarProtParameters(ϕ::Array{Float64,1},E::Float64)
    # Use asserts to ensure parameters a sensible
    @assert E >= 0.0 "Number of enzymes per cell cannot be negative"

    # Similar asserts for the vector parameters
    @assert length(ϕ) == 3 "Only considering 3 protein fractions"
    @assert all(ϕ .>= 0) "Protein fractions cannot be negative"
    @assert sum(ϕ) == 1 "Protein fractions must sum to one"
    return(VarProtParameters(ϕ,E))
end
