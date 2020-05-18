# Script that provides analytic functions such as Lyapunov exponents for use in the simulations
using Assembly

export Force, nForce

# function to return value of theta between 0 and 1
function bound_θ(S::Float64,P::Float64,T::Float64,η::Float64,ΔG0::Float64)
    # Catch perverse cases that sometimes arise
    if S <= 0.0
        θs = 1.0
    elseif P <= 0.0
        θs = 0.0
    else
        θs = Q(S,P)/Keq(T,η,ΔG0)
    end
    return(min(1,θs))
end

# function to return the rate in symbolic form
function symb_rate(S::Sym,θ::Sym)
    # Define kinetic/thermodynamic parameters
    qm, KS, kr, Kq = symbols("qm, KS, kr, Kq")
    # Start with numerator
    n = qm*S*(1 - θ)
    # Then find denominator
    d = KS + S*(1 + kr*θ)
    # Divide numerator by denominator
    q = n/d
    return(q)
end

# function to find the forces symbolically
function Force(ps::InhibParameters,F::Array{Sym,1})
    # Check Jacobian provided is the right size
    @assert length(F) == ps.N+ps.M "Preallocated force vector incorrect size"
    m, g, N, η, κ, M, δ = symbols("m, g, N, η, κ, M, δ")
    # And reaction level symbols
    qm, KS, kr, Kq = symbols("qm, KS, kr, Kq")
    # Loop over all microbes to begin with
    for i = 1:ps.N
        ex = -m
        # Loop over reactions
        for j = 1:ps.mics[i].R
            # Find the jth reaction
            rc = ps.reacs[ps.mics[i].Reacs[j]]
            # Make symbols for substrate and product
            S = symbols("M$(rc.Rct)")
            θ = symbols("$(rc.Rct)θ$(rc.Prd)N$i")
            # Then make symbolic rate
            q = symb_rate(S,θ)
            # Multiple expression by relevant factors
            ex += η*q
            # Now sub in all the reaction level parameters
            ex = subs(ex,qm=>ps.mics[i].qm[j],KS=>ps.mics[i].KS[j])
            ex = subs(ex,kr=>ps.mics[i].kr[j],η=>ps.mics[i].η[j])
        end
        ex *= g*N
        # Sub in strain/environment level variables
        F[i] = subs(ex,g=>ps.mics[i].g,m=>ps.mics[i].m,N=>"N$(i)")
    end
    # Then loop over metabolites
    for i = ps.N+1:ps.N+ps.M
        # Set metabolite number
        Mn = i-ps.N
        # Add external dynamics
        ex = κ - δ*M
        # Sub in relevant values
        ex = subs(ex,κ=>ps.κ[Mn],δ=>ps.δ[Mn],M=>"M$(Mn)")
        # Loop over all strains
        for j = 1:ps.N
            # Then loop over each reaction for each strain
            for k = 1:ps.mics[j].R
                # Find reactant metabolite number
                nR = ps.reacs[ps.mics[j].Reacs[k]].Rct
                # Find product metabolite number
                nP = ps.reacs[ps.mics[j].Reacs[k]].Prd
                # Check if metabolite is used as reactant
                if nR == Mn
                    # Make symbols for substrate and product
                    S = symbols("M$(nR)")
                    θ = symbols("$(nR)θ$(nP)N$j")
                    # Then make symbolic rate
                    q = symb_rate(S,θ)
                    # Then multiply by population
                    ex -= q*N
                    # sub in microbe identifier
                    ex = subs(ex,"N"=>"N$(j)")
                    # Then sub in all metabolite level parameters
                    ex = subs(ex,qm=>ps.mics[j].qm[k],KS=>ps.mics[j].KS[k])
                    ex = subs(ex,kr=>ps.mics[j].kr[k])
                # Or if metabolite is produced as product
                elseif nP == Mn
                    # Make symbols for substrate and product
                    S = symbols("M$(nR)")
                    θ = symbols("$(nR)θ$(nP)N$j")
                    # Then make symbolic rate
                    q = symb_rate(S,θ)
                    # Then multiply by population
                    ex += q*N
                    # sub in microbe identifier
                    ex = subs(ex,"N"=>"N$(j)")
                    # Then sub in all metabolite level parameters
                    ex = subs(ex,qm=>ps.mics[j].qm[k],KS=>ps.mics[j].KS[k])
                    ex = subs(ex,kr=>ps.mics[j].kr[k])
                end
            end
        end
        # No further substitutions required
        F[i] = ex
    end
    return(F)
end

# function to find numerical values of forces at a particular point in the space
function nForce(F::Array{Sym,1},C::Array{Float64,1},ps::InhibParameters)
    # Copy F to a new object to prevent overwriting
    f = copy(F)
    # Sub in steady state population values
    for i = 1:ps.N
        for j = 1:ps.N+ps.M
            f[j] = subs(f[j],"N$i"=>C[i])
        end
    end
    # Sub in steady state concentrations values
    for i = ps.N+1:ps.N+ps.M
        for j = 1:ps.N+ps.M
            f[j] = subs(f[j],"M$(i-ps.N)"=>C[i])
        end
    end
    # Finally sub in theta values
    for i = 1:ps.N
        for j = 1:ps.mics[i].R
            rc = ps.reacs[ps.mics[i].Reacs[j]]
            θ = bound_θ(C[ps.N+rc.Rct],C[ps.N+rc.Prd],ps.T,ps.mics[i].η[j],rc.ΔG0)
            for k = 1:ps.N+ps.M
                f[k] = subs(f[k],"$(rc.Rct)θ$(rc.Prd)N$(i)"=>θ)
            end
        end
    end
    # convert vector into a float
    f = convert(Array{Float64},f)
    return(f)
end
