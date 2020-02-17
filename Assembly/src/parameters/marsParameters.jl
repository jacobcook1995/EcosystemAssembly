# The function and struct used in this script were adapted from code orginally written by Tom Clegg

# I WILL NEED TO ALTER THE NAMES AND DEFINITIONS TO SUIT MY PURPOSES

"""
    Parameters(N::Int64,M::Int64,u::Array{Float64,2},m::Vector{Float64},g::Vector{Float64},l::Vector{Float64},h::Vector{Float64})
Type containing the parameters for a simuation.
# Arguments
- `N::Int64`: Number of Consumers
- `M::Int64`: Number of Resourcces
- `c::Array{Float64,2}`: Array of consumer preference of resource i by consumer j
- `m::Vector{Float64}`: Vector of maintenance energy costs
- `g::Vector{Float64}`: Vector of proportionality constants between growth and biomass
- `l::Vector{Float64}`: Vector of leakage fraction of particular resources
- `h::Vector{Float64}`: Vector of external resource supply rates
- `w::Vector{Float64}`: Vector of resource energy values
"""
struct Parameters
    N::Int64;
    M::Int64;
    c::Array{Float64,2}
    m::Vector{Float64}
    g::Vector{Float64}
    l::Vector{Float64}
    h::Vector{Float64}
    w::Vector{Float64}
end

"""
    make_Parameters(N::Int64,M::Int64,c::Array{Float64,2},m::Vector{Float64},g::Vector{Float64},l::Vector{Float64},h::Vector{Float64},w::Vector{Float64})
Helper function used internally. Takes values for parameters and returns a `Parameters`
object. Also does checks internally to make sure the values are correct.
"""
function make_Parameters(N::Int64,M::Int64,c::Array{Float64,2},m::Vector{Float64},g::Vector{Float64},l::Vector{Float64},h::Vector{Float64},w::Vector{Float64})

# Use asserts to ensure that arrays are correct sizes
@assert size(c) == (N,M) "Consumer preference array (c) is the wrong size"
@assert length(m) == N "Vector of maintenance costs (m) is the wrong length"
@assert length(g) == M "Vector of proportionality constants (g) vector is the wrong length"
@assert length(l) == M "Vector of leakage fractions (l) is wrong length"
@assert length(h) == M "Vector of external resource supplies (h) is the wrong length"
@assert length(w) == M "Vector of resource values (w) is the wrong length "

# Use asserts to avoid unphysical proportions
@assert all(l .<= 1) "One or more l values are > 1"

return(Parameters(N,M,c,m,g,l,h,w))
end
