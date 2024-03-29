# This script of the parameters for the full model, including proteome parameters
# This script depends on the Reaction type that is defined before

# export very simple functions
export ↦

# Overloadable version of getfields
↦(s, f) = getfield(s, f)

export Reaction, make_Reaction

export MicrobeP, make_MicrobeP, FullParameters, make_FullParameters, extinction

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
    ID::Int64
    Rct::Int64
    Prd::Int64
    ΔG0::Float64
end

"""
    make_Reaction(ID::Int64,Rct::Int64,Prd::Int64,ΔG0::Float64)
Helper function used internally. Takes values for parameters and returns a `Reaction`object.
Also does checks internally to make sure the values are correct.
"""
function make_Reaction(ID::Int64, Rct::Int64, Prd::Int64, ΔG0::Float64)
    @assert ID>0 "Reaction must be given a positive ID"
    @assert Rct>0 "All reactants must have positive IDs"
    @assert Prd>0 "All products must have positive IDs"
    @assert Prd!=Rct "Reactions cannot have same reactant and product"

    return (Reaction(ID, Rct, Prd, ΔG0))
end

"""
    MicrobeP(MC::Int64,γm::Float64,ρ::Float64,Kγ::Float64,Pb::Float64,d::Float64,ϕH::Float64,KΩ::Float64,
    fd::Float64,R::Int64,Reacs::Vector{Int64},η::Vector{Float64},kc::Vector{Float64},KS::Vector{Float64},
    kr::Vector{Float64},n::Vector{Int64},ϕP::Vector{Float64})
Type containing the parameters for a particular microbial strain, this now includes proteome parameters.
# Arguments
- `MC::Int64`: Mass of cell in amino acids.
- `γm::Float64`: Maximum elongation rate in amino acids per minute (per ribosome)
- `ρ::Float64`: ATP per sythesis step
- `Kγ::Float64`: Threshold energy in ATP per cell
- `Pb::Float64`: Proportion of ribosomes bound
- `d::Float64`: Biomass loss rate (~death rate)
- `ϕH::Float64`: Fraction of proteome allocated to housekeeping proteins
- `KΩ::Float64`: Saturation constant for proteome fraction with energy
- `fd::Float64`: Number of doublings required to reorder proteome
- `R::Int64`: Number of reactions
- `Reacs::Vector{Int64}`: Reaction numbers
- `η::Vector{Float64}`: ATP generated per mole of reaction
- `kc::Vector{Float64}`: Catalytic rate constants
- `KS::Vector{Float64}`: Substrate saturation constants
- `kr::Vector{Float64}`: reversibility factors
- `n::Vector{Int64}`: Array of number of amino acids per protein type, [r,p,h]
- `ϕP::Vector{Float64}`: Vector of metabolic fractions
"""
struct MicrobeP
    MC::Int64
    γm::Float64
    ρ::Float64
    Kγ::Float64
    Pb::Float64
    d::Float64
    ϕH::Float64
    KΩ::Float64
    fd::Float64
    R::Int64
    Reacs::Vector{Int64}
    η::Vector{Float64}
    kc::Vector{Float64}
    KS::Vector{Float64}
    kr::Vector{Float64}
    n::Vector{Int64}
    ϕP::Vector{Float64}
end

"""
    make_MicrobeP(MC::Int64,γm::Float64,ρ::Float64,Kγ::Float64,Pb::Float64,d::Float64,ϕH::Float64,KΩ::Float64,
    fd::Float64,R::Int64,Reacs::Vector{Int64},η::Vector{Float64},kc::Vector{Float64},KS::Vector{Float64},
    kr::Vector{Float64},n::Vector{Int64},ϕP::Vector{Float64})
Helper function used internally. Takes values for parameters and returns a `MicrobeP` object.
Also does checks internally to make sure the values are correct.
"""
function make_MicrobeP(MC::Int64, γm::Float64, ρ::Float64, Kγ::Float64, Pb::Float64,
                       d::Float64,
                       ϕH::Float64, KΩ::Float64, fd::Float64, R::Int64,
                       Reacs::Vector{Int64},
                       η::Vector{Float64}, kc::Vector{Float64}, KS::Vector{Float64},
                       kr::Vector{Float64},
                       n::Vector{Int64}, ϕP::Vector{Float64})
    # Check that physical parameters have been provided
    @assert R>0 "Number of reactions must be positive"
    @assert MC>0 "Cell mass must be positive"
    @assert γm>=0.0 "Elongation rate cannot be negative"
    @assert ρ>0.0 "ATP per synthesis step must be positive"
    @assert Kγ>0.0 "Threshold energy must be positive"
    @assert d>0.0 "Death rate must be positive"
    @assert 0.0<=Pb<=1.0 "Proportion of ribosomes bound has to be between 0 and 1"
    @assert 0.0<=ϕH<=1.0 "Housekeeping proteins have to be between 0 and 100% of total"
    @assert KΩ>0.0 "Proteome fraction saturation constant must be positive"
    @assert fd>0.0 "Number of doublings required must be positive"
    # Check that the vectors have the right length
    @assert length(Reacs)==R "Vector of reactions is the wrong length"
    @assert length(η)==R "Vector of ATP generation rates (η) is the wrong length"
    @assert length(kc)==R "Vector of maximal rates (qm) is the wrong length"
    @assert length(KS)==R "Vector of saturation constants (KS) is the wrong length"
    @assert length(kr)==R "Vector of reversibilities (kr) is the wrong length"
    @assert length(ϕP)==R "Vector of metabolic fractions is the wrong length"
    # Check that values in these vectors are plausible
    @assert all(η .>= 0.0) "η values cannot reasonably be negative"
    @assert all(kc .>= 0.0) "Catalytic rate constants cannot be negative"
    @assert all(KS .>= 0.0) "Saturation constants cannot be negative"
    @assert all(kr .>= 0.0) "Cannot have negative values for the reversibility"
    @assert all(Reacs .> 0) "All reactions are supposed to be indicated by a positive numbers"
    @assert length(n)==3 "Only considering 3 protein types"
    @assert all(n .> 0) "All proteins must have positive mass"
    @assert all(ϕP .>= 0.0) "All metabolic fractions must be non-negative"
    @assert sum(ϕP)≈1.0 "Metabolic fractions should sum to 1"
    return (MicrobeP(MC, γm, ρ, Kγ, Pb, d, ϕH, KΩ, fd, R, Reacs, η, kc, KS, kr, n, ϕP))
end

"""
    FullParameters(N::Int64,M::Int64,O::Int64,T::Float64,κ::Vector{Float64},δ::Vector{Float64},
    reacs::Vector{Reaction},mics::Vector{MicrobeP})
Type containing the parameters for a simuation.
# Arguments
- `N::Int64`: Number of Consumers
- `M::Int64`: Number of Resources
- `O::Int64`: Number of Reactions
- `T::Float64`: Temperature of system
- `κ::Vector{Float64}`: Vector of external resource supply rates
- `δ::Vector{Float64}`: Vector of resource decay rates
- `reacs::Vector{Reaction}`: Vector of reactions
- `mics::Vector{MicrobeP}`: Vector of microbial strains
"""
struct FullParameters
    N::Int64
    M::Int64
    O::Int64
    T::Float64
    κ::Vector{Float64}
    δ::Vector{Float64}
    reacs::Vector{Reaction}
    mics::Vector{MicrobeP}
end

"""
    make_FullParameters(N::Int64,M::Int64,O::Int64,T::Float64,κ::Vector{Float64},δ::Vector{Float64},
    reacs::Vector{Reaction},mics::Vector{MicrobeP})
Helper function used internally. Takes values for parameters and returns a `FullParameters`object.
Also does checks internally to make sure the values are correct.
"""
function make_FullParameters(N::Int64, M::Int64, O::Int64, T::Float64, κ::Vector{Float64},
                             δ::Vector{Float64},
                             reacs::Vector{Reaction}, mics::Vector{MicrobeP})
    # Use asserts to ensure that arrays are correct sizes
    @assert length(κ)==M "Vector of external resource supplies (κ) is the wrong length"
    @assert length(δ)==M "Vector of decay rates (δ) is the wrong length"
    @assert length(reacs)==O "Vector of reactions is the wrong length"
    @assert length(mics)==N "Vector of microbes is the wrong length"

    # Use asserts to avoid unphysical values
    @assert all(δ .>= 0) "One or more decay rates in δ are negative"
    @assert T>0.0 "Temperature must be positive"

    # Check that reaction numbering is reasonable
    @assert ((reacs .↦ :ID)==collect(1:O)) "Reaction IDs must be sequential from 1 upwards"
    @assert all((reacs .↦ :Rct) .<= M) "Reactant numbers cannot exceed number of metabolites"
    @assert all((reacs .↦ :Prd) .<= M) "Product numbers cannot exceed number of metabolites"
    # Check that no reaction is repeated
    for i in 1:O
        @assert all(((reacs .↦ :Rct)[i] .!= (reacs .↦ :Rct)[1:end .!= i]) .|
                    ((reacs .↦ :Prd)[i] .!= (reacs .↦ :Prd)[1:end .!= i])) "Reaction $i is repeated"
    end

    # Check that the reactions given in the vector of microbes exist in the vector of reactions
    for i in 1:N
        @assert all((mics .↦ :Reacs)[i] .<= O) "Microbe $i assigned to reaction that doesn't exist"
    end
    return (FullParameters(N, M, O, T, κ, δ, reacs, mics))
end

"""
    extinction(ps::FullParameters,ext::BitArray{1})
Helper function used internally. Recreate a FullParameters with extinct strains removed.
"""
function extinction(ps::FullParameters, ext::BitArray{1})
    # The majority of the parameter set is unaffected by the extinctions
    M = ps.M
    O = ps.O
    T = ps.T
    κ = ps.κ
    δ = ps.δ
    reacs = ps.reacs
    # Remove extinct strains from the count
    N = ps.N - sum(ext)
    # Preallocate array of new microbes
    mics = Array{MicrobeP, 1}(undef, N)
    # Now loop over microbes to retain the surviving strains
    k = 0
    for i in eachindex(ext)
        # Check if not marked as extinct
        if ext[i] != 1
            k += 1
            mics[k] = ps.mics[i]
        end
    end
    return (make_FullParameters(N, M, O, T, κ, δ, reacs, mics))
end
