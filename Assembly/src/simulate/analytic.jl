# Script that provides analytic functions such as Lyapunov exponents for use in the simulations
using Assembly

export Jacobian, Jacobian_test

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

# Function to make substrate differentiated function
function Sdiff()
    # Make symbols and standard function
    S, P = symbols("S, P")
    q = symb_rate(S,P)
    # Differentiate with respect to substrate
    dq = diff(q,S)
    return(dq)
end

# Function to make product differentiated function
function Pdiff()
    # Make symbols and standard function
    S, P = symbols("S, P")
    q = symb_rate(S,P)
    # Differentiate with respect to product
    dq = diff(q,P)
    return(dq)
end

# function to make the Jacobian for a particular parameter set
function Jacobian(ps::InhibParameters,J::Array{Sym,2})
    # Check Jacobian provided is the right size
    @assert size(J) == (ps.N+ps.M,ps.N+ps.M) "Preallocated Jacobian incorrect size"
    # Create symbols for repeated used variables
    m, g, d, η, N = symbols("m, g, d, η, N")
    # Define kinetic/thermodynamic parameters
    qm, KS, kr, Kq = symbols("qm, KS, kr, Kq")
    # Make differentiated expressions to sub into later
    diffS = Sdiff()
    diffP = Pdiff()
    for j = 1:size(J,2)
        for i = 1:size(J,1)
            # Deal with diagonals seperatly
            if i == j
                if i <= ps.N
                    # First add constant term
                    ex = -m*g
                    # Add all the J^grow values
                    for k = 1:ps.mics[i].R
                        # Find the kth reaction
                        rc = ps.reacs[ps.mics[i].Reacs[k]]
                        # Make symbols for substrate and product
                        S = symbols("M$(rc.Rct)")
                        P = symbols("M$(rc.Prd)")
                        # Then make symbolic rate
                        q = symb_rate(S,P)
                        # Multiple expression by relevant factors
                        ex += g*η*q
                        # Sub in all reaction level parameters
                        ex = subs(ex,qm=>ps.mics[i].qm[k],KS=>ps.mics[i].KS[k])
                        ex = subs(ex,kr=>ps.mics[i].kr[k],η=>ps.mics[i].η[k])
                        # Find value of equilbrium constant
                        K = Keq(ps.T,ps.mics[i].η[k],ps.reacs[ps.mics[i].Reacs[k]].ΔG0)
                        # And then sub it in
                        ex = subs(ex,Kq=>K)
                    end
                    # Now subsitute for all microbe level variables
                    J[i,j] = subs(ex,g=>ps.mics[i].g,m=>ps.mics[i].m)
                # Metabolite diagonals
                else
                    ex = -d
                    # Find metabolite number
                    Mn = j - ps.N
                    # Loop over all microbes
                    for k = 1:ps.N
                        # Then loop over each reaction for each microbe
                        for l = 1:ps.mics[k].R
                            # Find reactant metabolite number
                            nR = ps.reacs[ps.mics[k].Reacs[l]].Rct
                            # Find product metabolite number
                            nP = ps.reacs[ps.mics[k].Reacs[l]].Prd
                            # If metabolite is a reactant
                            if Mn == nR
                                ex -= diffS*N
                                # sub in metabolite identifiers
                                ex = subs(ex,"S"=>"M$(Mn)","P"=>"M$(nP)")
                                # Then sub in all metabolite level parameters
                                ex = subs(ex,qm=>ps.mics[k].qm[l],KS=>ps.mics[k].KS[l])
                                # Find value of equilbrium constant
                                K = Keq(ps.T,ps.mics[k].η[l],ps.reacs[ps.mics[k].Reacs[l]].ΔG0)
                                # Then sub it in
                                ex = subs(ex,kr=>ps.mics[k].kr[l],Kq=>K)
                            # If metabolite is a product
                            elseif Mn == nP
                                ex += diffP*N
                                # sub in metabolite identifiers
                                ex = subs(ex,"S"=>"M$(nR)","P"=>"M$(Mn)")
                                # Then sub in all metabolite level parameters
                                ex = subs(ex,qm=>ps.mics[k].qm[l],KS=>ps.mics[k].KS[l])
                                # Find value of equilbrium constant
                                K = Keq(ps.T,ps.mics[k].η[l],ps.reacs[ps.mics[k].Reacs[l]].ΔG0)
                                # Then sub it in
                                ex = subs(ex,kr=>ps.mics[k].kr[l],Kq=>K)
                            end
                        end
                        # Sub in specific population identifer
                        ex = subs(ex,N=>"N$k")
                    end
                    J[i,j] = subs(ex,d=>ps.δ[i-ps.N])
                end
            # All non diagonal pop cases are zero
            elseif i <= ps.N && j <= ps.N
                J[i,j] = 0
            # Then do non-diagonal metabolite cases
            elseif i > ps.N && j > ps.N
                ex = 0.0
                # Find metabolite number
                Mn = i - ps.N # metabolite changed
                Mnd = j - ps.N # metabolite to differentiate by
                # Loop over all microbes
                for k = 1:ps.N
                    # Then loop over each reaction for each microbe
                    for l = 1:ps.mics[k].R
                        # Find reactant metabolite number
                        nR = ps.reacs[ps.mics[k].Reacs[l]].Rct
                        # Find product metabolite number
                        nP = ps.reacs[ps.mics[k].Reacs[l]].Prd
                        # If metabolite is a reactant
                        if Mnd == nR && Mn == nP
                            ex += diffS*N
                            # sub in metabolite identifiers
                            ex = subs(ex,"S"=>"M$(Mnd)","P"=>"M$(Mn)")
                            # Then sub in all metabolite level parameters
                            ex = subs(ex,qm=>ps.mics[k].qm[l],KS=>ps.mics[k].KS[l])
                            # Find value of equilbrium constant
                            K = Keq(ps.T,ps.mics[k].η[l],ps.reacs[ps.mics[k].Reacs[l]].ΔG0)
                            # Then sub it in
                            ex = subs(ex,kr=>ps.mics[k].kr[l],Kq=>K)
                        # If metabolite is a product
                        elseif Mnd == nP && Mn == nR
                            ex -= diffP*N
                            # sub in metabolite identifiers
                            ex = subs(ex,"S"=>"M$(Mn)","P"=>"M$(Mnd)")
                            # Then sub in all metabolite level parameters
                            ex = subs(ex,qm=>ps.mics[k].qm[l],KS=>ps.mics[k].KS[l])
                            # Find value of equilbrium constant
                            K = Keq(ps.T,ps.mics[k].η[l],ps.reacs[ps.mics[k].Reacs[l]].ΔG0)
                            # Then sub it in
                            ex = subs(ex,kr=>ps.mics[k].kr[l],Kq=>K)
                        end
                    end
                    # Sub in specific population identifer
                    ex = subs(ex,N=>"N$k")
                end
                J[i,j] = ex
            # Population differentiated by metabolite case
            elseif i <= ps.N && j > ps.N
                ex = 0.0
                Mnd = j - ps.N
                # Sum over reactions
                for k = 1:ps.mics[i].R
                    # Find reactant metabolite number
                    nR = ps.reacs[ps.mics[i].Reacs[k]].Rct
                    # Find product metabolite number
                    nP = ps.reacs[ps.mics[i].Reacs[k]].Prd
                    # Check if differentiable by substrate
                    if nR == Mnd
                        ex += η*diffS
                        # sub in metabolite identifiers
                        ex = subs(ex,"S"=>"M$(Mnd)","P"=>"M$(nP)")
                        # Then sub in all metabolite level parameters
                        ex = subs(ex,qm=>ps.mics[i].qm[k],KS=>ps.mics[i].KS[k])
                        # Find value of equilbrium constant
                        K = Keq(ps.T,ps.mics[i].η[k],ps.reacs[ps.mics[i].Reacs[k]].ΔG0)
                        # Then sub it in
                        ex = subs(ex,kr=>ps.mics[i].kr[k],Kq=>K,η=>ps.mics[i].η[k])
                    # Or differntiable by product
                    elseif nP == Mnd
                        ex += η*diffP
                        # sub in metabolite identifiers
                        ex = subs(ex,"S"=>"M$(nR)","P"=>"M$(Mnd)")
                        # Then sub in all metabolite level parameters
                        ex = subs(ex,qm=>ps.mics[i].qm[k],KS=>ps.mics[i].KS[k])
                        # Find value of equilbrium constant
                        K = Keq(ps.T,ps.mics[i].η[k],ps.reacs[ps.mics[i].Reacs[k]].ΔG0)
                        # Then sub it in
                        ex = subs(ex,kr=>ps.mics[i].kr[k],Kq=>K,η=>ps.mics[i].η[k])
                    end
                end
                # Multiple expression by constant factors
                ex *= g*N
                J[i,j] = subs(ex,N=>"N$(i)",g=>ps.mics[i].g)
            # Metabolite differntiated by population case
            else
                ex = 0
                # Find meatbolite number
                Mn = i - ps.N
                # sum over all reactions for species j
                for k = 1:ps.mics[j].R
                    # Find reactant metabolite number
                    nR = ps.reacs[ps.mics[j].Reacs[k]].Rct
                    # Find product metabolite number
                    nP = ps.reacs[ps.mics[j].Reacs[k]].Prd
                    # Check if reaction involves metabolite as substrate
                    if nR == Mn
                        # Make symbols for substrate and product
                        S = symbols("M$(nR)")
                        P = symbols("M$(nP)")
                        # Then make symbolic rate
                        q = symb_rate(S,P)
                        ex -= q
                        # Sub in all reaction level parameters
                        ex = subs(ex,qm=>ps.mics[j].qm[k],KS=>ps.mics[j].KS[k])
                        ex = subs(ex,kr=>ps.mics[j].kr[k])
                        # Find value of equilbrium constant
                        K = Keq(ps.T,ps.mics[j].η[k],ps.reacs[ps.mics[j].Reacs[k]].ΔG0)
                        # And then sub it in
                        ex = subs(ex,Kq=>K)
                    # Or as product
                    elseif nP == Mn
                        # Make symbols for substrate and product
                        S = symbols("M$(nR)")
                        P = symbols("M$(nP)")
                        # Then make symbolic rate
                        q = symb_rate(S,P)
                        ex += q
                        # Sub in all reaction level parameters
                        ex = subs(ex,qm=>ps.mics[j].qm[k],KS=>ps.mics[j].KS[k])
                        ex = subs(ex,kr=>ps.mics[j].kr[k])
                        # Find value of equilbrium constant
                        K = Keq(ps.T,ps.mics[j].η[k],ps.reacs[ps.mics[j].Reacs[k]].ΔG0)
                        # And then sub it in
                        ex = subs(ex,Kq=>K)
                    end
                end
                J[i,j] = ex
            end
        end
    end
    return(J)
end

# Alternative Jacobian that makes a vector of dynamics and then directly differentiates
# This is easier to make correct but will be slower
# So I am adding this as a test
function Jacobian_test(ps::InhibParameters)
    # Make symbols that are going to be frequently used
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
    # Preallocate Jacobian matrix
    J = Array{Sym,2}(undef,ps.N+ps.M,ps.N+ps.M)
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
