# A script to numerically analyse a two step enzyme reaction
using Assembly
using SymPy
using Plots
import PyPlot

function twostep()
    # Define the symbols needed to define the process
    k1, k2, k3, k4, K1, K2, K3, K4 = symbols("k1, k2, k3, k4, K1, K2, K3, K4")
    A, B, C, E1, E10, E2, E20 = symbols("A, B, C, E1, E10, E2, E20")
    # Start with simple definitions of enzyme concentrations
    E1X = E10 - E1
    E2X = E20 - E2
    # Write out expressions for enzyme concentrations
    exr1 = E1*A*k1 + E1*B*K2 - k2*E1X - K1*E1X
    exr2 = E2*B*k3 + E2*C*K4 - k4*E2X - K3*E2X
    # Solve these expressions to find free enzyme concentrations at steady state
    E1 = solve(exr1,E1)
    E2 = solve(exr2,E2)
    # Want to remove the (only) elements from the arrays
    E1 = E1[1]
    E2 = E2[1]
    # Redefine the complexed enzymes
    E1X = E10 - E1
    E2X = E20 - E2
    # Now write out an expression for steady state B value
    exr3 = E1X*k2 + E2X*K3 - E2*B*k3 - E1*B*K2
    # Then solve for B
    Bs = solve(exr3,B)
    # Choose reasonable parameter values
    T = 293.15
    cE1 = 100000.0
    cE2 = cE1
    rk1 = 1000.0
    rk2 = 100.0
    rk3 = 1000.0
    rk4 = 100.0
    cA = 100.0
    cC = 1.0
    # Sub non-varying parameters into B
    Bs = subs.(Bs,E10=>cE1,E20=>cE2,A=>cA,C=>cC,k1=>rk1,k2=>rk2,k3=>rk3,k4=>rk4)
    # Overall equlbrium constant, and reaction quoitent
    Keq = 1.0
    Q = cC/cA
    # Make list of equilbrium constants to search
    Keq1 = [collect(0.0001:0.0001:0.001);collect(0.001:0.001:0.01);collect(0.01:0.01:0.1);collect(0.1:0.1:1.0);collect(1.0:1.0:10.0)]
    sec = [100.0,10.0,1.0,0.1,0.01]
    # Number of subdivisons to use
    N = length(Keq1)
    M = length(sec)
    # Preallocate vectors for plotting
    fs = zeros(N,M)
    G1 = zeros(N,M)
    # Loop over a range of equilbrium constant values
    for i = 1:N
        Keq2 = Keq/Keq1[i]
        for j = 1:M
            # Find combined reverse rate
            rA = (rk1*rk2)/(Keq1[i])
            rB = (rk3*rk4)/(Keq2)
            # Use vector of scenarios to choose relative reverse reaction rates
            rK1 = sqrt(rA/sec[j])
            rK2 = sec[j]*rK1
            # At the moment they correlate
            rK3 = sqrt(rB/sec[j])
            rK4 = sec[j]*rK3
            # Substitute them into an expression for B
            Bt = subs.(Bs,K1=>rK1,K2=>rK2,K3=>rK3,K4=>rK4)
            # Take the second +ve value
            cB = convert(Float64,Bt[2])
            # This expression is one way of calculating the reaction rate
            exr4 = E2*B*k3 - K3*E2X
            # Sub in all values to find value of the flux
            fs[i,j] = subs(exr4,B=>cB,C=>cC,E20=>cE2,k3=>rk3,K3=>rK3,k4=>rk4,K4=>rK4)
            # Now calculating reaction quotients for each step
            Q1 = cB/cA
            Q2 = cC/cB
            # Use to calculate free energy dissipated at each step
            G1[i,j] = Rgas*T*(log(Q1)-log(Keq1[i]))
        end
    end
    GT = Rgas*T*(log(Q)-log(Keq))
    # Now setup plotting
    pyplot()
    # Make labels
    lbs = Array{String,2}(undef,1,length(sec))
    for i = 1:length(sec)
        lbs[i] = "ratio = $(sec[i])"
    end
    # Set a color-blind friendly palette
    theme(:wong2,dpi=200)
    plot(G1/GT,fs,labels=lbs,ylabel="Overall reaction flux",xlabel="Fraction dissipated on first step")
    savefig("Output/50percent.png")
    return(nothing)
end

@time twostep()
