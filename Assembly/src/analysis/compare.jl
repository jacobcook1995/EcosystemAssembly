# A script to compare the output of the normal Marsland model and our alternative model
# with thermodynamic end product inhibition included
using Assembly
using Plots
import PyPlot
using JLD

# function to compare the two different η competition cases
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
    savefig("Output/Popsη1.png")
    plot(T1,C1[:,ps1.N+1:ps1.N+2],label="")
    savefig("Output/Concsη1.png")
    # Plot second case as well
    plot(T2,C2[:,1],label="consumer",color=ps2.N)
    for i = 2:ps2.N
        plot!(T2,C2[:,i],label="η = $(round(ps2.mics[i].η[1],digits=3))",color=i-1)
    end
    savefig("Output/Popsη2.png")
    plot(T2,C2[:,ps2.N+1:ps2.N+3],label="")
    savefig("Output/Concsη2.png")
    return(nothing)
end

# test function to run an example simulation
function compare()
    println("Successfully compiled.")
    # Want to investigate a chain of microbes
    N = 20
    M = 100
    O = 2*M
    mR = 5.0
    sdR = 1.0
    Tmax = 100.0
    mq = 1.0
    sdq = 0.1
    mK  = 0.1
    sdK = 0.01
    mk = 10.0
    sdk = 1.0
    for i = 1:50
        # Now make the parameter set
        ps = initialise(N,M,O,mR,sdR,mq,sdq,mK,sdK,mk,sdk)
        # save this parameter set
        jldopen("Temp/Paras/ps$(i).jld", "w") do file
            # addrequire(file, Assembly)
            write(file, "ps", ps)
        end
        @time C, T = inhib_simulate(ps,Tmax)
        # count number of species that have grown
        c = count(x->(x>=2.0), C[end,1:N])
        println("Simulation $i: $c species have grown")
        if i == 50
            # Run plotting
            pyplot(dpi=200)
            plot(T,C[:,1:N],label="")
            savefig("Output/TestPop.png")
            plot(T,C[:,N+1:end],label="")
            savefig("Output/TestConc.png")
            println(C[end,1:N])
        end
    end
    return(nothing)
end

# function to run test versions of scripts for the
function test()
    # Set simulation time
    Tmax = 25.0
    # Find filename to read in from argument
    ps = load("Temp/Paras/ps$(ARGS[1]).jld","ps")
    # remake parameter set with maintenance energy
    mics = Array{Microbe,1}(undef,ps.N)
    for i = 1:ps.N
        mics[i] = make_Microbe(1.0,ps.mics[i].g,ps.mics[i].R,ps.mics[i].Reacs,ps.mics[i].η,ps.mics[i].qm,ps.mics[i].KS,ps.mics[i].kr)
    end
    ps = make_InhibParameters(ps.N,ps.M,ps.O,ps.T,ps.κ,ps.δ,ps.reacs,mics)
    C, T = test_inhib_simulate(ps,Tmax)
    # If run is successful then do plotting
    pyplot(dpi=200)
    plot(T,C[:,1:ps.N],label="")
    savefig("Output/TestPop$(ARGS[1]).png")
    plot(T,C[:,ps.N+1:end],label="")
    savefig("Output/TestConc$(ARGS[1]).png")
    println(C[end,1:ps.N])
    return(nothing)
end

if length(ARGS) < 1
    @time compare()
else
    @time test()
end
