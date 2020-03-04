# Script that attempts to implement an alternative version of the Marsland model that
# utilizes themodynamic end product inhibition
using Assembly

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
