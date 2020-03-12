# The function and struct used in this script were adapted from code orginally written by Tom Clegg

# Export Parameters so that it can be used elsewhere in my code
export MarsParameters, make_MarsParameters, reduce_Parameters

"""
    MarsParameters(N::Int64,M::Int64,c::Array{Float64,2},m::Vector{Float64},g::Vector{Float64},l::Vector{Float64},
    κ::Vector{Float64},w::Vector{Float64},D::Array{Float64,2},δ::Vector{Float64})
Type containing the parameters for a simuation.
# Arguments
- `N::Int64`: Number of Consumers
- `M::Int64`: Number of Resourcces
- `c::Array{Float64,2}`: Array of consumer preference of resource i by consumer j
- `m::Vector{Float64}`: Vector of maintenance energy costs
- `g::Vector{Float64}`: Vector of proportionality constants between growth and biomass
- `l::Vector{Float64}`: Vector of leakage fraction of particular resources
- `κ::Vector{Float64}`: Vector of external resource supply rates
- `w::Vector{Float64}`: Vector of resource energy values
- `D::Array{Float64,2}`: Array of metabolic transistions
- `δ::Vector{Float64}`: Vector of resource decay rates
"""
struct MarsParameters
    N::Int64;
    M::Int64;
    c::Array{Float64,2}
    m::Vector{Float64}
    g::Vector{Float64}
    l::Vector{Float64}
    κ::Vector{Float64}
    w::Vector{Float64}
    D::Array{Float64,2}
    δ::Vector{Float64}
end

"""
    make_MarsParameters(N::Int64,M::Int64,c::Array{Float64,2},m::Vector{Float64},g::Vector{Float64},l::Vector{Float64},
    κ::Vector{Float64},w::Vector{Float64},D::Array{Float64,2},δ::Vector{Float64})
Helper function used internally. Takes values for parameters and returns a `Parameters`object.
Also does checks internally to make sure the values are correct.
"""
function make_MarsParameters(N::Int64,M::Int64,c::Array{Float64,2},m::Vector{Float64},g::Vector{Float64},l::Vector{Float64},κ::Vector{Float64},w::Vector{Float64},D::Array{Float64,2},δ::Vector{Float64})

# Use asserts to ensure that arrays are correct sizes
@assert size(c) == (N,M) "Consumer preference array (c) is the wrong size"
@assert length(m) == N "Vector of maintenance costs (m) is the wrong length"
@assert length(g) == N "Vector of proportionality constants (g) vector is the wrong length"
@assert length(l) == M "Vector of leakage fractions (l) is wrong length"
@assert length(κ) == M "Vector of external resource supplies (κ) is the wrong length"
@assert length(w) == M "Vector of resource values (w) is the wrong length "
@assert size(D) == (M,M) "Metabolic matrix (D) is the wrong size"
@assert length(δ) == M "Vector of decay rates (δ) is the wrong length"

# Use asserts to avoid unphysical values
@assert all(l .<= 1) "One or more l values are > 1"
@assert all(l .>= 0) "One or more l values are negative"
@assert all(sum(D,dims=1) .≈ 1) "One or more rows in D does not sum to 1"
@assert all(δ .>= 0) "One or more decay rates in δ are negative"

return(MarsParameters(N,M,c,m,g,l,κ,w,D,δ))
end

"""
    reduce_Parameters(ps::MarsParameters,rem::Array{Int64,1})
Helper function used internally. Takes an old parameter set and removes the elements specified in rem.
"""
function reduce_Parameters(ps::MarsParameters,rem::Array{Int64,1})

    # Use asserts to check that valid microbes for removal have been considered
    @assert length(unique(rem)) == length(rem) "Cannot remove a species twice"
    @assert all(0 .< rem .<= ps.N) "Cannot remove species that doesn't exist, species numbers run from 1 to $(ps.N)"

    # Determine number of species remaining
    nN = ps.N - length(rem)
    # Find indices of
    keep = collect(1:ps.N) .∉ rem
    nc = ps.c[keep,:]
    nm = ps.m[keep]
    ng = ps.g[keep]
    return(make_MarsParameters(nN,ps.M,nc,nm,ng,ps.l,ps.κ,ps.w,ps.D,ps.δ))
end