# A script to compare the output of the normal Marsland model and our alternative model
# with thermodynamic end product inhibition included
using Assembly
using Plots
import PyPlot

# function to compare Marsland and our extended models
function compare()
    println("Successfully compiled.")
    # Changing now to investigate η competition
    N = 11
    Tmax = 25.0
    mq = 1.0
    sdq = 0.1
    mK  = 0.1
    sdK = 0.01
    mk = 10.0
    sdk = 1.0
    ps = initialise_η2(N,mq,sdq,mK,sdK,mk,sdk)
    C, T = inhib_simulate(ps,Tmax)
    # Run plotting
    pyplot(dpi=200)
    plot(T,C[:,1],label="η = $(round(ps.mics[1].η[1],digits=3))") # Change this label at somepoint
    for i = 2:N
        plot!(T,C[:,i],label="η = $(round(ps.mics[i].η[1],digits=3))")
    end
    savefig("Output/PopTest.png")
    plot(T,C[:,N+1:N+3],label="")
    savefig("Output/ConcTest.png")
end

@time compare()
