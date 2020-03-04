# A script to compare the output of the normal Marsland model and our alternative model
# with thermodynamic end product inhibition included
using Assembly
using Plots
import PyPlot

# function to compare Marsland and our extended models
function compare()
    # Going to start with a small number of consumers and metabolities so that it runs fast, is easy to debug
    N = 20
    M = 100
    Tmax = 100.0
    O = 15
    C, T, ps = inhib_simulate(N,M,O,Tmax)
    # Run plotting
    pyplot(dpi=200)
    plot(T,C[:,1:N],label="")
    savefig("Output/PopTest.png")
    plot(T,C[:,N+1:N+M],label="")
    savefig("Output/ConcTest.png")
end

@time compare()
