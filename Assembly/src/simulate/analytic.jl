# Script that provides analytic functions such as Lyapunov exponents for use in the simulations
using Assembly

export Force, nForce, Jacobian, nJacobian

# function to return value of theta between 0 and 1
function bound_θ(S::Float64, P::Float64, T::Float64, η::Float64, ΔG0::Float64)
    # Catch perverse cases that sometimes arise
    if S <= 0.0
        θs = 1.0
    elseif P <= 0.0
        θs = 0.0
    else
        θs = Q(S, P) / Keq(T, η, ΔG0)
    end
    return (min(1, θs))
end

# function to return the rate in symbolic form
function symb_rate(S::SymPy.Sym, θ::SymPy.Sym)
    # Define kinetic/thermodynamic parameters
    kc, KS, kr = symbols("kc, KS, kr")
    # Start with numerator
    n = kc * S * (1 - θ)
    # Then find denominator
    d = KS + S * (1 + kr * θ)
    # Divide numerator by denominator
    q = n / d
    return (q)
end

# function to find the forces symbolically
function Force(ps::FullParameters, F::Array{SymPy.Sym, 1})
    # Check Force provided is the right size
    @assert length(F)==3 * ps.N + ps.M "Preallocated force vector incorrect size"
    # Define ecosystem and microbe level symbols
    N, η, κ, M, δ, d, γ, a, Kγ, Pb, ϕR, nr, np = symbols("N, η, κ, M, δ, d, γ, a, Kγ, fb, ϕR, nr, np")
    MC, ϕH, ρ, KΩ, fd = symbols("MC, ϕH, ρ, KΩ, fd")
    # And reaction level symbols
    kc, KS, kr, ϕP = symbols("kc, KS, kr, ϕP")
    # Loop over all microbes to begin with
    for i in 1:(ps.N)
        ex = -d
        ex += (γ * a / (a + Kγ)) * ϕR * Pb / nr
        ex *= N
        # Sub in strain/environment level variables
        ex = subs(ex, γ => ps.mics[i].γm, Kγ => ps.mics[i].Kγ, Pb => ps.mics[i].Pb,
                  nr => ps.mics[i].n[1])
        F[i] = subs(ex, d => ps.mics[i].d, a => "a$(i)", N => "N$(i)", ϕR => "ϕR$(i)")
    end
    # Then loop over metabolites
    for i in (ps.N + 1):(ps.N + ps.M)
        # Set metabolite number
        Mn = i - ps.N
        # Add external dynamics
        ex = κ - δ * M
        # Sub in relevant values
        ex = subs(ex, κ => ps.κ[Mn], δ => ps.δ[Mn], M => "M$(Mn)")
        # Loop over all strains
        for j in 1:(ps.N)
            # Then loop over each reaction for each strain
            for k in 1:(ps.mics[j].R)
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
                    q = symb_rate(S, θ)
                    # Then multiply by population, and amount of enzyme
                    ex -= q * N * MC * ϕP * (1 - ϕR - ϕH) / (np * NA)
                    # sub in microbe identifier and parameters
                    ex = subs(ex, N => "N$(j)", MC => ps.mics[j].MC, ϕR => "ϕR$(j)")
                    ex = subs(ex, ϕH => ps.mics[j].ϕH, np => ps.mics[j].n[2])
                    # Then sub in all metabolite level parameters
                    ex = subs(ex, kc => ps.mics[j].kc[k], KS => ps.mics[j].KS[k])
                    ex = subs(ex, kr => ps.mics[j].kr[k], ϕP => ps.mics[j].ϕP[k])
                    # Or if metabolite is produced as product
                elseif nP == Mn
                    # Make symbols for substrate and product
                    S = symbols("M$(nR)")
                    θ = symbols("$(nR)θ$(nP)N$j")
                    # Then make symbolic rate
                    q = symb_rate(S, θ)
                    # Then multiply by population, and amount of enzyme
                    ex += q * N * MC * ϕP * (1 - ϕR - ϕH) / (np * NA)
                    # sub in microbe identifier and parameters
                    ex = subs(ex, N => "N$(j)", MC => ps.mics[j].MC, ϕR => "ϕR$(j)")
                    ex = subs(ex, ϕH => ps.mics[j].ϕH, np => ps.mics[j].n[2])
                    # Then sub in all metabolite level parameters
                    ex = subs(ex, kc => ps.mics[j].kc[k], KS => ps.mics[j].KS[k])
                    ex = subs(ex, kr => ps.mics[j].kr[k], ϕP => ps.mics[j].ϕP[k])
                end
            end
        end
        # No further substitutions required
        F[i] = ex
    end
    # Next loop over the energy concentrations
    for i in (ps.N + ps.M + 1):(2 * ps.N + ps.M)
        # Find relevant microbe
        mic = ps.mics[i - ps.N - ps.M]
        # First set out growth rate
        ex = -(γ * a / (a + Kγ)) * Pb * ϕR / nr
        # Multiply it by translation cost and concentration
        ex *= (ρ * MC + a)
        # Then loop over reactions
        for j in 1:(mic.R)
            # Find the jth reaction
            rc = ps.reacs[mic.Reacs[j]]
            # Make symbols for substrate and product
            S = symbols("M$(rc.Rct)")
            θ = symbols("$(rc.Rct)θ$(rc.Prd)N$(i-ps.N-ps.M)")
            # Then make symbolic rate
            q = symb_rate(S, θ)
            # Multiple expression by relevant factors, including amount of enzyme
            ex += η * q * MC * ϕP * (1 - ϕR - ϕH) / (np)
            # Now sub in all the reaction level parameters
            ex = subs(ex, kc => mic.kc[j], KS => mic.KS[j], kr => mic.kr[j], η => mic.η[j],
                      ϕP => mic.ϕP[j])
        end
        # Sub in all the relevant variables
        ex = subs(ex, γ => mic.γm, a => "a$(i-ps.N-ps.M)", Kγ => mic.Kγ, Pb => mic.Pb,
                  nr => mic.n[1])
        F[i] = subs(ex, ϕR => "ϕR$(i-ps.N-ps.M)", ρ => mic.ρ, MC => mic.MC, ϕH => mic.ϕH,
                    np => mic.n[2])
    end
    # Finally do the same for the ribosome fractions
    for i in (2 * ps.N + ps.M + 1):(3 * ps.N + ps.M)
        # Find relevant microbe
        mic = ps.mics[i - 2 * ps.N - ps.M]
        # Subtract current ribosome fraction
        ex = -ϕR
        # Add the "optimal" ribosome fraction
        ex += (a / (a + KΩ)) * (1 - ϕH)
        # Then divide by the characteristic time scale
        ex /= fd * nr / ((γ * a / (a + Kγ)) * Pb * ϕR)
        # Sub in all the relevant variable
        ex = subs(ex, ϕR => "ϕR$(i-2*ps.N-ps.M)", a => "a$(i-2*ps.N-ps.M)", KΩ => mic.KΩ)
        F[i] = subs(ex, ϕH => mic.ϕH, fd => mic.fd, nr => mic.n[1], γ => mic.γm,
                    Kγ => mic.Kγ, Pb => mic.Pb)
    end
    return (F)
end

# function to find numerical values of forces at a particular point in the space
function nForce(F::Array{SymPy.Sym, 1}, C::Array{Float64, 1}, ps::FullParameters)
    # Copy F to a new object to prevent overwriting
    f = copy(F)
    # Sub in steady state population values
    for i in 1:(ps.N)
        for j in 1:(3 * ps.N + ps.M)
            f[j] = subs(f[j], "N$i" => C[i])
        end
    end
    # Sub in steady state concentrations values
    for i in (ps.N + 1):(ps.N + ps.M)
        for j in 1:(3 * ps.N + ps.M)
            f[j] = subs(f[j], "M$(i-ps.N)" => C[i])
        end
    end
    # Then sub in theta values
    for i in 1:(ps.N)
        for j in 1:(ps.mics[i].R)
            rc = ps.reacs[ps.mics[i].Reacs[j]]
            θ = bound_θ(C[ps.N + rc.Rct], C[ps.N + rc.Prd], ps.T, ps.mics[i].η[j], rc.ΔG0)
            for k in 1:(3 * ps.N + ps.M)
                f[k] = subs(f[k], "$(rc.Rct)θ$(rc.Prd)N$(i)" => θ)
            end
        end
    end
    # Next sub in energies
    for i in (ps.N + ps.M + 1):(2 * ps.N + ps.M)
        for j in 1:(3 * ps.N + ps.M)
            f[j] = subs(f[j], "a$(i-ps.N-ps.M)" => C[i])
        end
    end
    # Finally sub in ribosome fractions
    for i in (2 * ps.N + ps.M + 1):(3 * ps.N + ps.M)
        for j in 1:(3 * ps.N + ps.M)
            f[j] = subs(f[j], "ϕR$(i-2*ps.N-ps.M)" => C[i])
        end
    end
    # convert vector into a float
    f = convert(Array{Float64}, f)
    return (f)
end

# function to find an analytic form of the Jacobian
function Jacobian(F::Array{SymPy.Sym, 1}, J::Array{SymPy.Sym, 2}, ps::FullParameters,
                  cons::Array{Float64, 1})
    # Check Force provided is the right size
    @assert size(J)==(length(F), length(F)) "preallocated Jacobian wrong size"
    @assert size(J)==(3 * ps.N + ps.M, 3 * ps.N + ps.M) "preallocated Jacobian wrong size"
    @assert length(cons)==ps.M "concentrations wrong length"

    # Copy F to a new object to prevent overwriting
    f = copy(F)
    # Loop over strains to sub in theta values
    for i in 1:(ps.N)
        # For each microbe loop over the reactions it has
        for j in 1:(ps.mics[i].R)
            # Find reaction
            rc = ps.reacs[ps.mics[i].Reacs[j]]
            # And corresponding θ symbol
            θ = "$(rc.Rct)θ$(rc.Prd)N$(i)"
            # Find S and P
            S, P = symbols("M$(rc.Rct), M$(rc.Prd)")
            # Find equilibrium constant
            K = Keq(ps.T, ps.mics[i].η[j], rc.ΔG0)
            # Check if substrate concentration is zero
            if cons[rc.Rct] == 0.0
                # In this case substitute 1 in to remove the relevant terms
                for k in 1:(3 * ps.N + ps.M)
                    f[k] = subs(f[k], θ => 1.0)
                end
            else
                # Sub true theta value into the vector of forces
                for k in 1:(3 * ps.N + ps.M)
                    f[k] = subs(f[k], θ => P / (S * K))
                end
            end
        end
    end
    # Now loop over variables to differentiate by
    for j in 1:(3 * ps.N + ps.M)
        # Find variable to differentiate by
        if j <= ps.N
            dx = symbols("N$j")
        elseif j <= ps.N + ps.M
            dx = symbols("M$(j-ps.N)")
        elseif j <= 2 * ps.N + ps.M
            dx = symbols("a$(j-ps.N-ps.M)")
        else
            dx = symbols("ϕR$(j-2*ps.N-ps.M)")
        end
        for i in 1:(3 * ps.N + ps.M)
            J[i, j] = diff(f[i], dx)
        end
    end
    return (J)
end

# function to find numerical values of forces at a particular point in the space
function nJacobian(J::Array{SymPy.Sym, 2}, C::Array{Float64, 1}, ps::FullParameters)
    # Copy F to a new object to prevent overwriting
    Jc = copy(J)
    # Sub in steady state population values
    for i in 1:(ps.N)
        for j in 1:(3 * ps.N + ps.M)
            for k in 1:(3 * ps.N + ps.M)
                Jc[j, k] = subs(Jc[j, k], "N$i" => C[i])
            end
        end
    end
    # Sub in steady state concentrations values
    for i in (ps.N + 1):(ps.N + ps.M)
        for j in 1:(3 * ps.N + ps.M)
            for k in 1:(3 * ps.N + ps.M)
                Jc[j, k] = subs(Jc[j, k], "M$(i-ps.N)" => C[i])
            end
        end
    end
    # Next sub in energies
    for i in (ps.N + ps.M + 1):(2 * ps.N + ps.M)
        for j in 1:(3 * ps.N + ps.M)
            for k in 1:(3 * ps.N + ps.M)
                Jc[j, k] = subs(Jc[j, k], "a$(i-ps.N-ps.M)" => C[i])
            end
        end
    end
    # Finally sub in ribosome fractions
    for i in (2 * ps.N + ps.M + 1):(3 * ps.N + ps.M)
        for j in 1:(3 * ps.N + ps.M)
            for k in 1:(3 * ps.N + ps.M)
                Jc[j, k] = subs(Jc[j, k], "ϕR$(i-2*ps.N-ps.M)" => C[i])
            end
        end
    end
    # convert vector into a float
    Jc = convert(Array{Float64}, Jc)
    return (Jc)
end
