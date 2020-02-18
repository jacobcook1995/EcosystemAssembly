# Script that attempts to implement the model laid out in Marsland et al, PLoS Comput Biol 15: e1006793.
# This provides a point of comparison for our extended model
using Assembly
using Distributions

# function to find rate of intake of a particular resource, this is currently a linear function
function vin(pref::Float64,conc::Float64)
    return(pref*conc)
end

# function to find rate of output of a particular resource, this matches vin as a linear function
# wb is the value of the output resource, vins is 1D as only need inputs for specific microbe type
# Da is the vector of transformations from resource a to other resources
function vout(l::Array{Float64},wb::Float64,w::Array{Float64,1},vins::Array{Float64,1},Da::Array{Float64,1})
    # Initialise vout to zero
    vo = 0
    # Loop over all resourses
    for i = 1:length(l)
        vo += (w[i]/wb)*Da[i]*l[i]*vins[i]
    end
    return(vo)
end

# function to find the energy obtained for growth from all resources for a particular cell
function Jgrow(M::Int64,l::Array{Float64,1},vals::Array{Float64,1},vins::Array{Float64,1})
    Jg = 0
    # Loop over all resources
    for i = 1:M
        Jg += (1-l[i])*vals[i]*vins[i]
    end
    return(Jg)
end

# function to implement the consumer resource dynamics
# I currently provide vin and vout in order to avoid uneeded allocations => Is this sensible?
# Also think that pop and concs need to be merged => Do this later if neccesary
function dynamics(ps::Parameters,dxdt::Array{Float64,1},pop::Array{Float64,1},concs::Array{Float64,1},vins::Array{Float64,2},vouts::Array{Float64,2})
    # First find and store intake rates
    for j = 1:ps.M
        for i = 1:ps.N
            vins[i,j] = vin(ps.c[i,j],concs[j])
        end
    end
    # Then use to find output rates
    for j = 1:ps.M
        for i = 1:ps.N
            vouts[i,j] = vout(ps.l,ps.w[j],ps.w,vins[i,:],ps.D[j,:])
        end
    end
    # First consumer dynamics
    for i = 1:ps.N
        dxdt[i] = ps.g[i]*pop[i]*(Jgrow(ps.M,ps.l,ps.w,vins[i,:]) - ps.m[i])
    end
    # Then resource dynamics
    for i = ps.N+1:ps.N+ps.M
        # fist add external supply of resource and decay
        dxdt[i] = ps.κ[i-ps.N] - ps.δ[i-ps.N]
        # Then consider contribution by various microbes
        for j = 1:ps.N
            dxdt[i] += pop[i]*(vouts[j,i-ps.N] - vins[j,i-ps.N])
        end
    end
    # Finally return the new dxdt values
    return(dxdt)
end

# function to construct vector of metabolite types, very simple at momet but can tweak it if I wish
function Mtypes(M::Int64,Nt::Int64)
    # Simple error message here to avoid negative occurances of metabolite types
    @assert round(M/4) + Nt <= M "Cannot have more metabolite types than metabolites"
    # Initialise the vectors
    Ms = zeros(Int64,M)
    cM = zeros(Int64,1+Nt)
    # One in four metabolities are of prefered byproduct type
    cM[1] = round(M/4)
    # Otherwise divided evenly between byproduct types
    for i = 2:length(cM)-1
        cM[i] = round(3*M/(4*Nt))
    end
    # This adjusts for the rounding => Final catagory can often be bigger
    cM[end] = M - sum(cM[1:end-1])
    # Now need to do metabolite numbering
    ccM = accumulate(+,cM) # Find cumulative sum of cM
    count = 1
    while count <= M
        Ms[count] = findfirst(count.<=ccM)
        count += 1
    end
    return(Ms,cM)
end

# function to construct the metabolic matrix
function Dmatrix(M::Int64,fc::Float64,fs::Float64,d0::Float64,Ms::Array{Int64,1},cM::Array{Int64})
    # Checks that the entered parameters are reasonable
    @assert fc + fs <= 1 "Fractions cannot sum to greater than 1"
    @assert fc > 0 && fs > 0 "Fractions cannot be less than zero"
    @assert 0.0 < d0 <= 1.0 "Stochasticity parameter d0 must be greater than 0 and less than 1"
    # initialise matrix
    D = zeros(M,M)
    # preallocate vector for the parameters for the Dirichlet distribution
    Dps = zeros(M)
    # Sample each column from distribution
    for j = 1:M
        for i = 1:M
            if Ms[i] == 1 && Ms[j] == 1
                Dps[i] = d0*(fc+fs)/(cM[1])
            elseif Ms[i] == 1
                Dps[i] = d0*(fc)/(cM[1])
            elseif Ms[j] == 1
                Dps[i] = d0*(1-fc-fs)/(M-cM[1])
            elseif Ms[i] == Ms[j]
                Dps[i] = d0*fs/cM[Ms[j]]
            else
                Dps[i] = d0*(1-fs-fc)/(M-cM[Ms[j]]-cM[1])
            end
        end
        # Now need to use this column to sample from the Dirichlet distribution
        d = Dirichlet(Dps)
        D[:,j] .= rand(d)
    end
    return(D)
end

# function to run simulation of the Marsland model
function simulate()
    # Going to start with a small number of consumers and metabolities so that it runs fast, is easy to debug
    N = 5
    M = 20
    # All metabolities have the same value for simplicity
    w = ones(M)
    # And all proportionality constants are the same for simplicity
    g = ones(N)
    # All but resource 1 is not supplied
    κ = zeros(M)
    κ[1] = 100.0
    # Leakage fractions are assumed to be equal across all metabolities
    li = 0.1
    l = li*ones(M)
    # Assume that all δ's are equal
    δi = 1.0
    δ = δi*ones(M)
    # Find M types so that I can define c and D
    Nt = 3 # Number of additional types
    Ms, cM = Mtypes(M,Nt)
    # Find D using a function that samples from the Dirichlet distribution
    fc = 0.3 # These fractions are parameters that could be changed
    fs = 0.3
    d0 = 0.5 # Stochasticity parameter, worth fiddling with
    D = Dmatrix(M,fc,fs,d0,Ms,cM)
    # c is a bit more fiddly => Needs a function
    # m is fixed with a Guassian offset => Also needs to be a function
    return(nothing)
end

@time simulate()
