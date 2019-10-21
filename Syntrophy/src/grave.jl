# Graveyard for old functions that I should delete eventually

function gluc()
    # Nutrient variables
    α = 3.00*10^(-8)
    δ = 1.00*10^(-6) # Death rate that doesn't wash out cells
    # make nutrients
    # 1 = glucose, 2 = oxegen, 3 = bicarbonate, 4 = hydrogen ion
    nuts = [Nut(1,α,δ),Nut(2,6*α,δ),Nut(3,0,δ),Nut(4,0,δ)]
    # Now make reactions
    ΔG0 = -2843800.0
    reac = [React(1,[1,2,3,4],[-1,-6,6,6],ΔG0)]
    # microbe variables
    η1 = 41.2
    η2 = 41.0
    # maintainance equal
    m = 2.16*10^(-19)
    # both have same substrate and end product
    r = 1
    # make microbes
    mics = [Microbe(η1,m,r,δ),Microbe(η2,m,r,δ)]
    # Set intial populations and nutrient concentrations
    pops = 10.0*ones(length(mics))
    concs = zeros(length(nuts))
    concs[1] = 0.3# start with high amount of glucose
    concs[2] = 0.0018 # WHY NOT JUST FIX O2 and pH?
    concs[3] = 1.00*10^(-9)
    concs[4] = 1.00*10^(-7) # pH 7
    # Define some constants
    Y = 2.36*10^(13) # biomass yield cells per mol ATP
    Γ = 1.16*10^(12) # starvation rate cells per mol ATP (deficit)
    K = 2.00*10^(-10) # Saturation constant mol ml^(−1)
    qm = 4.44*10^(-13) # Maximal possible growth rate mol cell^(-1) s^(-1)
    ΔGATP = 75000.0 # Gibbs free energy of formation of ATP in a standard cell
    Temp = 312.0 # Temperature that growth is occuring at in Kelvin
    p = [Y,Γ,K,qm,ΔGATP,Temp]
    u0 = [concs;pops]
    tspan = (0.0,2500000.0)
    ex = [2,4]
    # Make reduced version of function inputting unchanging microbes
    f(du,u,p,t) = npops(du,u,p,nuts,reac,mics,ex,t)
    prob = ODEProblem(f,u0,tspan,p)
    sol = solve(prob,adaptive=false,dt=2000) # turned dt down to make plots look nicer
    # Now do plotting
    plot(sol.t,sol'[:,5:6],label=["\\eta = $((mics.↦:η)[1])","\\eta = $((mics.↦:η)[2])"])
    savefig("Output/Populations$(η1)vs$(η2).png")
    return(nothing)
end

# function that approximates solution as 1st order polynomial and then
function anstead1(sm::Float64,sH::Float64,sα::Float64,sδ::Float64,sη::Float64,sqm::Float64,skr::Float64,sKS::Float64,sKeq::Float64,sO2::Float64)
    P, S, X, m, η, δ, α, q = symbols("P S X m eta delta alpha q",positive=true)
    # Rate of change of substrate
    dS = -m*X/η - δ*S + α
    # Rate of change of product
    dP = 6*m*X/η - δ*P
    # Solve first expression for S
    X1 = SymPy.solve(dS,X)
    # Only one solution so overwrite
    X1 = X1[1]
    # Substitute into product rate of change to eliminate X
    dP = subs(dP,X=>X1)
    # Then try to solve this expression for P
    P1 = SymPy.solve(dP,P)
    # Only one solution so overwrite
    P1 = P1[1]
    # Now write out expression for rate of change of X
    dX = η*q - m
    # Need a lot more variables to define q
    qm, R, θ, KS, Q, Keq, kr, O2, H = symbols("qm R theta KS Q Keq kr O2 H",positive=true)
    # Now define q and sub expression in dX
    q1 = qm*R*(1-θ)/(KS+R*(1+kr*θ))
    dX = subs(dX,q=>q1)
    # Next sub for R
    R1 = S*O2^6
    dX = subs(dX,R=>R1)
    # Then need to define and expand θ and sub in to dX
    θ1 = Q/Keq
    Q1 = (P^6)*(H^6)/(S*(O2^6))
    θ1 = subs(θ1,Q=>Q1)
    dX = subs(dX,θ=>θ1)
    # Now sub in expression for to obtain expression in terms of S
    dX = subs(dX,P=>P1)
    # Complex expression, first simplify
    dX = simplify(dX)
    # Then remove denominator as solving equal to zero
    dX = dX*denom(dX)
    # Convert expression to a (6th order) polynomial in S
    dX = sympy.Poly(dX,S) # Needs to be lower case as this is qualifying the containing (python) module
    # Now reduce to first order polynomial and solve
    c0 = dX.coeffs()[7]
    c1 = dX.coeffs()[6]
    dX1 = c1*S + c0
    S1 = SymPy.solve(dX1,S)
    S1 = S1[1] # Take single solution
    # Then subsitute for the relevant values
    S1 = subs(S1,m=>sm,H=>sH,α=>sα,δ=>sδ,η=>sη,qm=>sqm,kr=>skr,KS=>sKS,Keq=>sKeq,O2=>sO2)
    # And do same for P
    P1 = subs(P1,S=>S1,α=>sα,δ=>sδ)
    # Now find X
    X1 = subs(X1,S=>S1,α=>sα,δ=>sδ,η=>sη,m=>sm)
    return(S1,P1,X1)
end

# overload function to look at how obtain more extensive (3rd order) approximate solution
function anstead3(sm::BigFloat,sH::BigFloat,sα::BigFloat,sδ::BigFloat,sη::BigFloat,sqm::BigFloat,skr::BigFloat,sKS::BigFloat,sKeq::BigFloat,sO2::BigFloat)
    P, S, X, m, η, δ, α, q = symbols("P S X m eta delta alpha q",positive=true)
    # Rate of change of substrate
    dS = -m*X/η - δ*S + α
    # Rate of change of product
    dP = 6*m*X/η - δ*P
    # Solve first expression for S
    X1 = SymPy.solve(dS,X)
    # Only one solution so overwrite
    X1 = X1[1]
    # Substitute into product rate of change to eliminate X
    dP = subs(dP,X=>X1)
    # Then try to solve this expression for P
    P1 = SymPy.solve(dP,P)
    # Only one solution so overwrite
    P1 = P1[1]
    # Now write out expression for rate of change of X
    dX = η*q - m
    # Need a lot more variables to define q
    qm, R, θ, KS, Q, Keq, kr, O2, H = symbols("qm R theta KS Q Keq kr O2 H",positive=true)
    # Now define q and sub expression in dX
    q1 = qm*R*(1-θ)/(KS+R*(1+kr*θ))
    dX = subs(dX,q=>q1)
    # Next sub for R
    R1 = S*O2^6
    dX = subs(dX,R=>R1)
    # Then need to define and expand θ and sub in to dX
    θ1 = Q/Keq
    Q1 = (P^6)*(H^6)/(S*(O2^6))
    θ1 = subs(θ1,Q=>Q1)
    dX = subs(dX,θ=>θ1)
    # Now sub in expression for to obtain expression in terms of S
    dX = subs(dX,P=>P1)
    # Complex expression, first simplify
    dX = simplify(dX)
    # Then remove denominator as solving equal to zero
    dX = dX*denom(dX)
    # Convert expression to a (6th order) polynomial in S
    dX = sympy.Poly(dX,S) # Needs to be lower case as this is qualifying the containing (python) module
    # Now reduce to third order polynomial and solve
    c0 = dX.coeffs()[7]
    c1 = dX.coeffs()[6]
    c2 = dX.coeffs()[5]
    c3 = dX.coeffs()[4]
    dX1 = c3*S^3 + c2*S^2 + c1*S + c0
    S1 = SymPy.solve(dX1,S)
    # Three solutions now, only first a real solution
    SA = subs(S1[1],m=>sm,H=>sH,α=>sα,δ=>sδ,η=>sη,qm=>sqm,kr=>skr,KS=>sKS,Keq=>sKeq,O2=>sO2)
    # And do same for P
    PA = subs(P1,S=>SA,α=>sα,δ=>sδ)
    # Now find X
    XA = subs(X1,S=>SA,α=>sα,δ=>sδ,η=>sη,m=>sm)
    # Convert S, P, X back to Float64's
    SA = convert(Float64,SA)
    PA = convert(Float64,PA)
    XA = convert(Float64,XA)
    return(SA,PA,XA)
end

# overload function to look at how obtain more extensive (4th order) approximate solution
function anstead4(sm::BigFloat,sH::BigFloat,sα::BigFloat,sδ::BigFloat,sη::BigFloat,sqm::BigFloat,skr::BigFloat,sKS::BigFloat,sKeq::BigFloat,sO2::BigFloat)
    P, S, X, m, η, δ, α, q = symbols("P S X m eta delta alpha q",positive=true)
    # Rate of change of substrate
    dS = -m*X/η - δ*S + α
    # Rate of change of product
    dP = 6*m*X/η - δ*P
    # Solve first expression for S
    X1 = SymPy.solve(dS,X)
    # Only one solution so overwrite
    X1 = X1[1]
    # Substitute into product rate of change to eliminate X
    dP = subs(dP,X=>X1)
    # Then try to solve this expression for P
    P1 = SymPy.solve(dP,P)
    # Only one solution so overwrite
    P1 = P1[1]
    # Now write out expression for rate of change of X
    dX = η*q - m
    # Need a lot more variables to define q
    qm, R, θ, KS, Q, Keq, kr, O2, H = symbols("qm R theta KS Q Keq kr O2 H",positive=true)
    # Now define q and sub expression in dX
    q1 = qm*R*(1-θ)/(KS+R*(1+kr*θ))
    dX = subs(dX,q=>q1)
    # Next sub for R
    R1 = S*O2^6
    dX = subs(dX,R=>R1)
    # Then need to define and expand θ and sub in to dX
    θ1 = Q/Keq
    Q1 = (P^6)*(H^6)/(S*(O2^6))
    θ1 = subs(θ1,Q=>Q1)
    dX = subs(dX,θ=>θ1)
    # Now sub in expression for to obtain expression in terms of S
    dX = subs(dX,P=>P1)
    # Complex expression, first simplify
    dX = simplify(dX)
    # Then remove denominator as solving equal to zero
    dX = dX*denom(dX)
    # Convert expression to a (6th order) polynomial in S
    dX = sympy.Poly(dX,S) # Needs to be lower case as this is qualifying the containing (python) module
    # Now reduce to first order polynomial and solve
    c0 = dX.coeffs()[7]
    c1 = dX.coeffs()[6]
    c2 = dX.coeffs()[5]
    c3 = dX.coeffs()[4]
    c4 = dX.coeffs()[3]
    dX1 = c4*S^4 + c3*S^3 + c2*S^2 + c1*S + c0
    S1 = SymPy.solve(dX1,S)
    # Four solutions now, first two have imaginary parts, fourth gives X and P negative
    SC = subs(S1[3],m=>sm,H=>sH,α=>sα,δ=>sδ,η=>sη,qm=>sqm,kr=>skr,KS=>sKS,Keq=>sKeq,O2=>sO2)
    # And do same for P
    PC = subs(P1,S=>SC,α=>sα,δ=>sδ)
    # Now find X
    XC = subs(X1,S=>SC,α=>sα,δ=>sδ,η=>sη,m=>sm)
    # Convert S, P, X back to Float64's
    SC = convert(Float64,SC)
    PC = convert(Float64,PC)
    XC = convert(Float64,XC)
    return(SC,PC,XC)
end

# This only works as an approximation whilst approximation to zero still works
# Function to find an approximate value for
function apθ(α::Float64,δ::Float64,S::Float64,H::Float64,O2::Float64,Keq::Float64)
    # Calculate ratio of constant products vs reactants
    rH = (H/O2)^6
    # Ratio of variable products vs reactants
    rS = ((6*((α/δ)-S))^6)/S
    # Use to calculate conventional form of θ
    θ = rH*rS/Keq
    # Print out values of all terms
    T1 = 6^6*(α/δ)^6
    θ = rH*T1/(Keq*S)
    return(θ)
end
