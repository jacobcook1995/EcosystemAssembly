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
    E0 = 100000.0
    rk1 = 1000.0
    rk2 = 100.0
    rk3 = 1000.0
    rk4 = 100.0
    rK1 = 10.0
    rK3 = 10.0
    cA = 100.0
    cC = 1.0
    # Sub non-varying parameters into B
    Bs = subs.(Bs,E10=>E0,E20=>E0,A=>cA,C=>cC,k1=>rk1,k2=>rk2,k3=>rk3,k4=>rk4,K1=>rK1,K3=>rK3)
    # Overall equlbrium constant
    Keq = 10.0
    # Preallocate vectors for plotting
    fs = zeros(1000)
    Ks = zeros(1000)
    # Loop over a range of equilbrium constant values
    for i = 1:1000
        Keq1 = i/10.0
        Keq2 = Keq/Keq1
        rK2 = (rk1*rk2)/(Keq1*rK1)
        rK4 = (rk3*rk4)/(Keq2*rK3)
        # Substitute them into an expression for B
        Bt = subs.(Bs,K2=>rK2,K4=>rK4)
        # Take the second +ve value
        cB = convert(Float64,Bt[2])
        # This expression is one way of calculating the reaction rate
        exr4 = E2*B*k3 - K3*E2X
        # Sub in all values to find value of the flux
        fs[i] = subs(exr4,B=>cB,C=>cC,E20=>E0,k3=>rk3,K3=>rK3,k4=>rk4,K4=>rK4)
        Ks[i] = Keq1
    end
    # Now setup plotting
    pyplot()
    # Set a color-blind friendly palette
    theme(:wong2,dpi=200)
    plot(Ks,fs)
    savefig("Output/test.png")
    return(nothing)
end

@time twostep()
