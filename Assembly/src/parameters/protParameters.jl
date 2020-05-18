# This script of the parameters for the single population proteome model
# This should be later integrated with inhibParameters.

export ProtParameters, make_ProtParameters

"""
    ProtParameters(MC::Int64,γm::Float64,T::Float64,n::Array{Int64,1})
Type containing the parameters for a simulation of the single population proteome model.
# Arguments
- `MC::Int64`: Mass of cell in amino acids.
- `γm::Float64`: Maximum elongation rate in amino acids per minute (per ribosome)
- `T::Float64`: Temperature that ecosystem is at
- `n::Array{Int64,1}`: Array of number of amino acids per protein type, [r,p,h]
"""
struct ProtParameters
    MC::Int64;
    γm::Float64;
    T::Float64;
    n::Array{Int64,1}
end

"""
    make_ProtParameters(MC::Int64,γm::Float64,T::Float64,n::Array{Int64,1})
Helper function used internally. Takes values for parameters and returns a `InhibParameters`object.
Also does checks internally to make sure the values are correct.
"""
function make_ProtParameters(MC::Int64,γm::Float64,T::Float64,n::Array{Int64,1})
    # Use asserts to ensure that arrays are correct sizes
    @assert MC > 0 "Cell mass must be positive"
    @assert γm >= 0.0 "Elongation rate cannot be negative"
    @assert T > 0.0 "Temperature must be positive"

    # Similar asserts for the vector parameters
    @assert length(n) == 3 "Only considering 3 protein types"
    @assert all(n .> 0) "All proteins must have positive mass"

    return(ProtParameters(MC,γm,T,n))
end
