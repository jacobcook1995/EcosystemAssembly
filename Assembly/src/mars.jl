# Script that attempts to implement the model laid out in Marsland et al, PLoS Comput Biol 15: e1006793.
# This provides a point of comparison for our extended model
using Assembly
using Distributions
using DifferentialEquations
using BenchmarkTools # REMOVE THIS ONCE I'VE DONE TESTING

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

# function to find the specialism of the consumers
function special(N::Int64,Nt::Int64)
    # Initialise the vectors
    Mp = zeros(Int64,N)
    # One more entry than number of metabolities types
    # [1] => generalist, [2] => specialist for metabolite etc
    cMp = zeros(Int64,1+Nt)
    # Half the species are generalists
    cMp[1] = round(N/2)
    # Randomly assign the remaining consumers generalisms
    for i = cMp[1]+1:N
        Mp[i] = rand(1:Nt)
    end
    # Then sort so that it is easier to read
    Mp = sort(Mp)
    # Count these and add them to the vector
    for i = 1:Nt
        cMp[i+1] = count(x->x==i,Mp)
    end
    return(Mp,cMp)
end

# function to construct the matrix of consumer preferences
function cmatrix(N::Int64,M::Int64,Mp::Array{Int64,1},cM::Array{Int64,1},Ms::Array{Int64,1},c0::Float64,c1::Float64,μc::Float64,qA::Float64)
    @assert qA >= 0.0 "Consumer preference strength parameter cannot be negative"
    @assert c0 >= 0.0 && c1 >= 0.0 "Cannot have negative preference coefficients"
    # Initialise matrix
    c = zeros(N,M)
    # Determine low value and high values
    l = c0/M
    h = c0/M + c1
    # Loop over all metabolites and all consumers
    for j = 1:M
        for i = 1:N
            # Check if species is a generalist
            if Mp[i] == 0
                p = μc/(M*c1)
            # Or a specialist-specialism pair
            elseif Mp[i] == Ms[j]
                p = (μc/(M*c1))*(1+qA*((M-cM[Mp[i]])/cM[Mp[i]]))
            else
                p = μc/(M*c1)*(1-qA)
            end
            # Draw random number if less than p then set peference high
            r = rand()
            if r <= p
                c[i,j] = h
            else
                c[i,j] = l
            end
        end
    end
    return(c)
end

# function to generate a vector of values for the maintenance energy requirments m
function mvector(N::Int64,mm::Float64,sdm::Float64)
    @assert mm - 5*sdm >= 0.0 "This choice could result in negative energy requirements"
    # Initialise vector of m
    m = zeros(N)
    # Make required Gaussian distribution using the provided mean (mm) and SD (sdm)
    d = Normal(mm,sdm)
    for i = 1:N
        m[i] = rand(d)
    end
    return(m)
end

# function to implement the consumer resource dynamics
# I currently provide vin and vout in order to avoid uneeded allocations => Is this sensible?
# ARE STATIC ARRAYS A GOOD IDEA? => Definetly are if I am not adding or removing species
function dynamics!(dx::Array{Float64,1},x::Array{Float64,1},t::Float64,ps::Parameters,vins::Array{Float64,2},vouts::Array{Float64,2})
    # First find and store intake rates
    for j = 1:ps.M
        for i = 1:ps.N
            vins[i,j] = vin(ps.c[i,j],x[ps.N+j])
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
        dx[i] = ps.g[i]*x[i]*(Jgrow(ps.M,ps.l,ps.w,vins[i,:]) - ps.m[i])
    end
    # Then resource dynamics
    for i = ps.N+1:ps.N+ps.M
        # fist add external supply of resource and decay
        dx[i] = ps.κ[i-ps.N] - ps.δ[i-ps.N]
        # Then consider contribution by various microbes
        for j = 1:ps.N
            dx[i] += x[i]*(vouts[j,i-ps.N] - vins[j,i-ps.N])
        end
    end
    # Finally return the new dxdt values
    return(dxdt)
end

# function to run simulation of the Marsland model
function simulate()
    # Going to start with a small number of consumers and metabolities so that it runs fast, is easy to debug
    N = 6
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
    Nt = 4 # Number of metabolite types
    Ms, cM = Mtypes(M,Nt)
    # Find D using a function that samples from the Dirichlet distribution
    fc = 0.3 # These fractions are parameters that could be changed
    fs = 0.3
    d0 = 0.2 # Stochasticity parameter, worth fiddling with
    D = Dmatrix(M,fc,fs,d0,Ms,cM)
    # Find generalism or specialism of consumers so that c can be found
    Mp, cMp = special(N,Nt)
    # Specify high and low expression levels
    hi = 1.00
    lo = 0.01
    # Average of 10% of metabolites selected
    μc = 0.1*M
    # Define c0 and c1 such that they result in correct expression values
    c0 = 0.01
    c1 = 1.00
    qA = 0.5 # preference strength parameter, worth fiddling with
    # Find c using a function that samples from a binary probability distribution
    c = cmatrix(N,M,Mp,cM,Ms,c0,c1,μc,qA)
    # Find m using a function that gives a Guassian offset
    mm = 1.0
    sdm = 0.1
    m = mvector(N,mm,sdm)
    # Now make the parameter set
    ps = make_Parameters(N,M,c,m,g,l,κ,w,D,δ)
    # Initialise vectors of concentrations and populations
    pop = ones(N)
    conc = zeros(M)
    x0 = [pop;conc]
    vins = zeros(N,M)
    vouts = zeros(N,M)
    # Give the function the parameters + preallocated memory
    dyns!(dx,x,t) = dynamics!(dx,x,t,ps,vins,vouts)
    # Try this now
    tspan = (0.0,100.0)
    prob = ODEProblem(dyns!,x0,tspan)
    println(prob)
    @benchmark solve(prob,Tsit5())
    return(nothing)
end

@time simulate()
