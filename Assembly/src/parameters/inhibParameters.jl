# This script is similar to MarsParameters but is for the inhibitory model

export InhibParameters, convert_Parameters

"""
    InhibParameters(N::Int64,M::Int64)
Type containing the parameters for a simuation.
# Arguments
- `N::Int64`: Number of Consumers
- `M::Int64`: Number of Resources
- `c::Array{Float64,2}`: Array of consumer preference of resource i by consumer j
- `m::Vector{Float64}`: Vector of maintenance energy costs
- `g::Vector{Float64}`: Vector of proportionality constants between growth and biomass
- `κ::Vector{Float64}`: Vector of external resource supply rates
- `δ::Vector{Float64}`: Vector of resource decay rates
"""
struct InhibParameters
    N::Int64;
    M::Int64;
    c::Array{Float64,2}
    m::Vector{Float64}
    g::Vector{Float64}
    κ::Vector{Float64}
    δ::Vector{Float64}
end

"""
   convert_Parameters(ps::MarsParameters)
Helper function used internally. Takes a MarsParameters object and converts it to an equivalent InhibParameters object.
"""
function convert_Parameters(ps::MarsParameters)
    # N and M unchanged
    N = ps.N
    M = ps.M
    # The consumer preference array should not change
    c = ps.c
    # Neither should the vectors g, m, κ and δ
    m = ps.m
    g = ps.g
    κ = ps.κ
    δ = ps.δ
    return(make_InhibParameters(N,M,c,m,g,κ,δ))
end

"""
    make_InhibParameters(N::Int64,M::Int64,c::Array{Float64,2},m::Vector{Float64},g::Vector{Float64},κ::Vector{Float64},δ::Vector{Float64})
Helper function used internally. Takes values for parameters and returns a `Parameters`object.
Also does checks internally to make sure the values are correct.
"""
function make_InhibParameters(N::Int64,M::Int64,c::Array{Float64,2},m::Vector{Float64},g::Vector{Float64},κ::Vector{Float64},δ::Vector{Float64})

# Use asserts to ensure that arrays are correct sizes
@assert size(c) == (N,M) "Consumer preference array (c) is the wrong size"
@assert length(m) == N "Vector of maintenance costs (m) is the wrong length"
@assert length(g) == N "Vector of proportionality constants (g) vector is the wrong length"
@assert length(κ) == M "Vector of external resource supplies (κ) is the wrong length"
@assert length(δ) == M "Vector of decay rates (δ) is the wrong length"

# Use asserts to avoid unphysical values
@assert all(δ .>= 0) "One or more decay rates in δ are negative"

return(InhibParameters(N,M,c,m,g,κ,δ))
end
