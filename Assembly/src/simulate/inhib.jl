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
    return(q)
end

# function to find the energy obtained for growth from all resources for a particular cell
function Jgrow()

    return(Jg)
end

# function to implement the consumer resource dynamics
# NEEDS TO BE SUBSTANTIALLY ALTERED
function dynamics!(dx::Array{Float64,1},x::Array{Float64,1},ps::MarsParameters,vins::Array{Float64,2},vouts::Array{Float64,2},t::Float64)
    # First find and store intake rates
    for j = 1:ps.M
        for i = 1:ps.N
            # Set value of vin to zero if it has gone negative
            vins[i,j] = max(vin(ps.c[i,j],x[ps.N+j]),0.0)
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
        # FINE APART FROM Jgrow WHICH NEEDS TO BE CHANGED
        dx[i] = ps.g[i]*x[i]*(Jgrow() - ps.m[i])
    end
    # Then resource dynamics
    for i = ps.N+1:ps.N+ps.M
        # fist add external supply of resource and decay
        dx[i] = ps.κ[i-ps.N] - ps.δ[i-ps.N]*x[i]
        # Then consider contribution by various microbes
        for j = 1:ps.N
            dx[i] += x[j]*(vouts[j,i-ps.N] - vins[j,i-ps.N])
        end
    end
    # Finally return the new dxdt values
    return(dx)
end
