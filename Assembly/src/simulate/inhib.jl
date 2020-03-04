# Script that attempts to implement an alternative version of the Marsland model that
# utilizes themodynamic end product inhibition
using Assembly

export inhib_simulate

# function to run simulation of the Marsland model
function initialise(N::Int64,M::Int64,O::Int64)
    # Assume that temperature T is constant at 20°C
    T = 293.15
    # And all proportionality constants are the same for simplicity
    g = ones(N)
    # All but resource 1 is not supplied
    κ = zeros(M)
    κ[1] = 100.0
    # Assume that all δ's are equal
    δi = 1.0
    δ = δi*ones(M)
    # Find m using a function that gives a Guassian offset
    mm = 1.0
    sdm = 0.1
    m = mvector(N,mm,sdm)
    # Preallocate vector of recations
    reacs = Array{Reaction,1}(undef,O)
    for i = 1:O
        reacs[i] = make_Reaction(i,Rct::Int64,Prd::Int64,ΔG0::Float64)
    end
    # Preallocate vector of microbes
    mics = Array{Microbe,1}(undef,N)
    # Then construct microbes
    for i = 1:N
        mics[i] = make_Microbe(m[i],g[i],R::Int64,Reacs::Vector{Int64},η::Vector{Float64},qm::Vector{Float64},KS::Vector{Float64},kr::Vector{Float64})
    end
    # Now make the parameter set
    ps = make_InhibParameters(N,M,O,T,κ,δ,reacs,mics)
    return(ps)
end

# function to find the reaction quotient Q, in the case of 1 to 1 stochiometery
function Q(S::Float64,P::Float64)
    Q = P/S
    return(Q)
end

# function to find the equilbrium constant
function Keq(T::Float64,η::Float64,ΔG0::Float64)
    Keq = exp((-ΔG0-η*ΔGATP)/(Rgas*T))
    return(Keq)
end

# function to find the thermodynamic term θ, for the case of 1 to 1 stochiometry
function θ(S::Float64,P::Float64,T::Float64,η::Float64,ΔG0::Float64)
    θ = Q(S,P)/Keq(T,η,ΔG0)
    return(θ)
end

# function to find the rate of substrate consumption by a particular reaction
function qs(S::Float64,P::Float64,T::Float64,η::Float64,ΔG0::Float64,qm::Float64,KS::Float64,kr::Float64)
    θs = θ(S,P,T,η,ΔG0)
    q = qm*S*(1-θs)/(KS + S*(1+kr*θs))
    # Ensure that negative value cannot be returned
    return(max(q,0.0))
end


# function to implement the consumer resource dynamics
function dynamics!(dx::Array{Float64,1},x::Array{Float64,1},ps::InhibParameters,rate::Array{Float64,2},t::Float64)
    # loop over the reactions to find reaction rate for each reaction for each strain
    for j = 1:ps.O
        # Find substrate and product for this reaction
        S = x[ps.N+ps.reacs[j].Rct]
        P = x[ps.N+ps.reacs[j].Prd]
        for i = 1:ps.N
            # Check if microbe i performs reaction j
            if j ∈ ps.mic[i].Reacs
                # Find index of this reaction in microbe
                k = findfirst(x->x==j,ps.mic[i].Reacs)
                # Use k to select correct kinetic parameters
                rate[i,j] = qs(S,P,ps.T,ps.mics[i].η[k],ps.reacs[j].ΔG0,ps.mics[i].qm[k],ps.mics[i].KS[k],ps.mics[i].kr[k])
            else
                rate[i,j] = 0.0
            end
        end
    end
    # Now want to use the rate matrix in the consumer dynamics
    for i = 1:ps.N
        # subtract maintenance
        dx[i] = -ps.mics.m[i]
        # Add all the non zero contributions
        for j = 1:ps.mics.R
            dx[i] += ps.mics.η[j]*rate[i,ps.mics.Reacs[j]]
        end
        # multiply by population and proportionality constant
        dx[i] *= ps.mics[i].g*x[i]
    end
    # Do basic resource dynamics
    for i = ps.N+1:ps.N+ps.M
        # fist add external supply of resource and decay
        dx[i] = ps.κ[i-ps.N] - ps.δ[i-ps.N]*x[i]
    end
    # Then loop over microbes
    for i = 1:ps.N
        # Loop over reactions for specific microbe
        for j = 1:ps.mics.R
            # Increase the product
            dx[ps.N+ps.reacs[ps.mics.Reacs[j]].Prd] += rate[i,ps.mics.Reacs[j]]*x[i]
            # and decrease the reactant
            dx[ps.N+ps.reacs[ps.mics.Reacs[j]].Rct] -= rate[i,ps.mics.Reacs[j]]*x[i]
        end
    end
    return(dx)
end

# Simulation code to run one instatnce of the simulation
# N is number of microbial strains, M is the number of metabolites
# O is the number of reactions, Tmax is how long the simulation is run for
function inhib_simulate(N::Int64,M::Int64,O::Int64,Tmax::Float64)
    # Make random parameter set of this size
    ps = initialise(N,M,O)
    # Initialise vectors of concentrations and populations
    pop = ones(N)
    conc = 0.1*ones(M) # Initial trace amount of each metabolite
    # Preallocate memory
    rate = zeros(N,O)
    # Now substitute preallocated memory in
    dyns!(dx,x,ps,t) = dynamics!(dx,x,ps,rate,t)
    # Find time span for this step
    tspan = (0,Tmax)
    x0 = [pop;conc]
    # Then setup and solve the problem
    prob = ODEProblem(dyns!,x0,tspan,ps)
    sol = solve(prob,isoutofdomain=(y,p,t)->any(x->x<0,y))
    return(sol',sol.t,ps)
end
