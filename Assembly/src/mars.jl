# Script that attempts to implement the model laid out in Marsland et al, PLoS Comput Biol 15: e1006793.
# This provides a point of comparison for our extended model
using Assembly

# function to find rate of intake of a particular resource, this is currently a linear function
function vin(pref::Float64,conc::Float64)
    return(pref*conc)
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
function dynamics(ps::Parameters,dxdt::Array{Float64,1},pop::Array{Float64,1},concs::Array{Float64,1},vins::Array{Float64,2},vouts::Array{Float64,2})
    # First find and store intake rates
    for j = 1:ps.M
        for i = 1:ps.N
            vins[i,j] = vin(ps.c[i,j],concs[j])
        end
    end
    # Then use to find output rates, NEED TO PUT A FORMULA IN HERE

    # First consumer dynamics
    for i = 1:ps.N
        dxdt[i] = ps.g[i]*pop[i]*(Jgrow(ps.M,ps.l,ps.w,vins[i,:]) - ps.m[i])
    end
    # Then resource dynamics
    for i = ps.N+1:ps.N+ps.M
        # fist add external supply of resource
        dxdt[i] = ps.h[i-ps.N]
        # Then consider contribution by various microbes
        for j = 1:ps.N
            dxdt[i] += pop[i]*(vouts[j,i] - vins[j,i])
        end
    end
    # Finally return the new dxdt values
    return(dxdt)
end
