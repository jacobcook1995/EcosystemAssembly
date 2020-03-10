# A script to compare the output of the normal Marsland model and our alternative model
# with thermodynamic end product inhibition included
using Assembly
using Plots
import PyPlot

# function to compare Marsland and our extended models
function ηcomparison()
    println("Successfully compiled.")
    # Changing now to investigate η competition
    N = 10
    Tmax = 100.0
    mq = 1.0
    mK  = 0.1
    mk = 10.0
    # Run alternative simulatiosn of the two cases
    ps1 = initialise_η(N,mq,mK,mk)
    C1, T1 = inhib_simulate(ps1,Tmax)
    # One extra strain as consumer exists
    ps2 = initialise_η2(N+1,mq,mK,mk)
    C2, T2 = inhib_simulate(ps2,Tmax)
    # Run plotting
    pyplot(dpi=200)
    plot(T1,C1[:,1],label="η = $(round(ps1.mics[1].η[1],digits=3))") # WRONG LABEL
    for i = 2:ps1.N
        plot!(T1,C1[:,i],label="η = $(round(ps1.mics[i].η[1],digits=3))")
    end
    savefig("Output/PopTest1.png")
    plot(T1,C1[:,ps1.N+1:ps1.N+2],label="")
    savefig("Output/ConcTest1.png")
    # Plot second case as well
    plot(T2,C2[:,1],label="consumer",color=ps2.N)
    for i = 2:ps2.N
        plot!(T2,C2[:,i],label="η = $(round(ps2.mics[i].η[1],digits=3))",color=i-1)
    end
    savefig("Output/PopTest2.png")
    plot(T2,C2[:,ps2.N+1:ps2.N+3],label="")
    savefig("Output/ConcTest2.png")
end

@time ηcomparison()
