# Script that provides analytic functions such as Lyapunov exponents for use in the simulations
using Assembly

export Jacobian, Lyapunov

# function to return the rate in symbolic form
function symb_rate(S::Sym,P::Sym)
    # Define kinetic/thermodynamic parameters
    qm, KS, kr, Kq = symbols("qm, KS, kr, Kq")
    # Start with numerator
    n = qm*S*(1 - P/(S*Kq))
    # Then find denominator
    d = KS + S*(1+(kr*P)/(S*Kq))
    # Divide numerator by denominator
    q = n/d
    return(q)
end

# function to make the Jacobian for a particular parameter set
function Jacobian(ps::InhibParameters,J::Array{Sym,2})
    # Check Jacobian provided is the right size
    @assert size(J) == (ps.N+ps.M,ps.N+ps.M) "Preallocated Jacobian incorrect size"
    m, g, N, η, κ, M, δ = symbols("m, g, N, η, κ, M, δ")
    # And reaction level symbols
    qm, KS, kr, Kq = symbols("qm, KS, kr, Kq")
    # Make empty vector
    f = Array{Sym,1}(undef,ps.N+ps.M)
    for i = 1:ps.N
        ex = -m
        # Loop over reactions
        for j = 1:ps.mics[i].R
            # Find the jth reaction
            rc = ps.reacs[ps.mics[i].Reacs[j]]
            # Make symbols for substrate and product
            S = symbols("M$(rc.Rct)")
            P = symbols("M$(rc.Prd)")
            # Then make symbolic rate
            q = symb_rate(S,P)
            # Multiple expression by relevant factors
            ex += η*q
            # Now sub in all the reaction level parameters
            ex = subs(ex,qm=>ps.mics[i].qm[j],KS=>ps.mics[i].KS[j])
            ex = subs(ex,kr=>ps.mics[i].kr[j],η=>ps.mics[i].η[j])
            # Find value of equilbrium constant
            K = Keq(ps.T,ps.mics[i].η[j],ps.reacs[ps.mics[i].Reacs[j]].ΔG0)
            # And then sub it in
            ex = subs(ex,Kq=>K)
        end
        ex *= g*N
        # Sub in strain/envioment level variables
        f[i] = subs(ex,g=>ps.mics[i].g,m=>ps.mics[i].m,N=>"N$(i)")
    end
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
                    P = symbols("M$(nP)")
                    # Then make symbolic rate
                    q = symb_rate(S,P)
                    # Then multiply by population
                    ex -= q*N
                    #  sub in metabolite and microbe identifiers
                    ex = subs(ex,"S"=>"M$(Mn)","P"=>"M$(nP)","N"=>"N$(j)")
                    # Then sub in all metabolite level parameters
                    ex = subs(ex,qm=>ps.mics[j].qm[k],KS=>ps.mics[j].KS[k])
                    # Find value of equilbrium constant
                    K = Keq(ps.T,ps.mics[j].η[k],ps.reacs[ps.mics[j].Reacs[k]].ΔG0)
                    # Then sub it in
                    ex = subs(ex,kr=>ps.mics[j].kr[k],Kq=>K)
                # Or if metabolite is produced as product
                elseif nP == Mn
                    # Make symbols for substrate and product
                    S = symbols("M$(nR)")
                    P = symbols("M$(nP)")
                    # Then make symbolic rate
                    q = symb_rate(S,P)
                    # Then multiply by population
                    ex += q*N
                    # sub in metabolite and microbe identifiers
                    ex = subs(ex,"S"=>"M$(nR)","P"=>"M$(Mn)","N"=>"N$(j)")
                    # Then sub in all metabolite level parameters
                    ex = subs(ex,qm=>ps.mics[j].qm[k],KS=>ps.mics[j].KS[k])
                    # Find value of equilbrium constant
                    K = Keq(ps.T,ps.mics[j].η[k],ps.reacs[ps.mics[j].Reacs[k]].ΔG0)
                    # Then sub it in
                    ex = subs(ex,kr=>ps.mics[j].kr[k],Kq=>K)
                end
            end
        end
        # No further substitutions required
        f[i] = ex
    end
    # Loop over variables to differntiate by
    for j = 1:ps.N+ps.M
        # Find variable to differentiate by
        if j <= ps.N
            dx = symbols("N$j")
        else
            dx = symbols("M$(j-ps.N)")
        end
        for i = 1:ps.N+ps.M
            J[i,j] = diff(f[i],dx)
        end
    end
    return(J)
end

# function to find the local eigenvalues and vectors for a partcilar Jacobian
function Lyapunov(J::Array{Sym,2},C::Array{Float64,1},ps::InhibParameters)
    # Check that Jacobian does correspond to
    @assert length(C) == size(J,1) "Steady state vector wrong size for Jacobian"
    for i = ps.N+1:ps.N+ps.M
        # Overwrite zero concentrations
        if C[i] == 0.0
            C[i] = 1e-10
        end
    end
    # Sub in steady state population values
    for i = 1:ps.N
        for j = 1:ps.N+ps.M
            for k = 1:ps.N+ps.M
                J[j,k] = subs(J[j,k],"N$i"=>C[i])
            end
        end
    end
    # Sub in steady state concentrations values
    for i = ps.N+1:ps.N+ps.M
        for j = 1:ps.N+ps.M
            for k = 1:ps.N+ps.M
                J[j,k] = subs(J[j,k],"M$(i-ps.N)"=>C[i])
            end
        end
    end
    # convert fully substituted Jacobian into a float
    J = convert(Array{Float64},J)
    # Then find eigenvalues
    λ = eigvals(J)
    v = eigvecs(J)
    return(λ,v)
end
