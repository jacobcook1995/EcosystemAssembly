# THIS IS A TEST SCRIPT DELETE WHEN DONE!!!
# DELETE EVENTUALLY
using Distributions

# function to construct vector of metabolite types, very simple at momet but can tweak it if I wish
function Mtypes(M::Int64,Nt::Int64)
    # Simple error message here to avoid negative occurances of metabolite types
    @assert round(M/4) + (Nt-1) <= M "Cannot have more metabolite types than metabolites"
    # Initialise the vectors
    Ms = zeros(Int64,M)
    cM = zeros(Int64,Nt)
    # One in four metabolities are of prefered byproduct type
    cM[1] = round(M/4)
    # Otherwise divided evenly between byproduct types
    for i = 2:length(cM)-1
        cM[i] = round(3*M/(4*(Nt-1)))
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

function test()
    # Choose parameters
    M = 20
    Nt = 4
    fc = 0.3 # These fractions are parameters that could be changed
    fs = 0.3
    d0 = 0.3 # Stochasticity parameter, worth fiddling with
    # Make M types
    Ms, cM = Mtypes(M,Nt)
    # That make stochiometric matrix
    D = Dmatrix(M,fc,fs,d0,Ms::Array{Int64,1},cM::Array{Int64})
    return(nothing)
end

@time test()
