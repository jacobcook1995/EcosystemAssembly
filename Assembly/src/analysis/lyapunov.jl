# A script to run Lyapunov analysis of our inhibition model.
# KEEP AS AN EXAMPLE SCRIPT FOR NOW
using Assembly
using JLD
using SymPy
using Plots
import PyPlot

# function to set up and run a lyapunov analysis of a simple case
function test_run()
    # only a small number of strains and metabolites for this analysis
    N = 4
    M = 4
    O = 4
    mR = 2.0
    sdR = 0.0
    mq = 1.0
    sdq = 0.1
    mK  = 0.1
    sdK = 0.01
    mk = 10.0
    sdk = 1.0
    Tmax = 100.0
    # Now make the parameter set, can just use the normal function for now
    ps = initialise_test(mq,mK,mk)
    # save this parameter set
    jldopen("Temp/Paras/psL.jld","w") do file
        write(file,"ps",ps)
    end
    # Preallocate Jacobian
    J = Array{Sym,2}(undef,ps.N+ps.M,ps.N+ps.M)
    # Find Jacobian using function
    J = Jacobian(ps,J)
    # Now run simulation
    C, T = inhib_simulate(ps,Tmax)
    # If run is successful then do plotting
    pyplot(dpi=200)
    plot(T,C[:,1:ps.N],label="")
    savefig("Output/LyaPop$(ARGS[1]).png")
    plot(T,C[:,ps.N+1:end],label="")
    savefig("Output/LyaConc$(ARGS[1]).png")
    return(nothing)
end

# function to run test versions of scripts for the
function lya()
    # Set simulation time
    Tmax = 100.0
    # Find filename to read in from argument
    ps = load("Temp/Paras/psL.jld","ps")
    # Preallocate Jacobian
    J = Array{Sym,2}(undef,ps.N+ps.M,ps.N+ps.M)
    # Find Jacobian using function
    J = Jacobian(ps,J)
    # Then run test simulation of it
    C, T = inhib_simulate(ps,Tmax)
    # Use final simulation results to find local Lyapunov exponents
    λ, v = Lyapunov(J,C[end,:],ps)
    for i = 1:length(λ)
        println("λ$i = $(λ[i])")
        println("v$i = $(v[:,i])")
    end
    # Add microbe and test if this changes things
    mic = make_Microbe(1.0,1.0,1,[3],[3.5],[1.0],[0.1],[10.0])
    ps = add_Microbe(ps,mic)
    # Preallocate Jacobian
    J = Array{Sym,2}(undef,ps.N+ps.M,ps.N+ps.M)
    # Find Jacobian using function
    J = Jacobian(ps,J)
    λ, v = Lyapunov(J,[C[end,1:ps.N-1];0.0;C[end,ps.N:end]],ps)
    for i = 1:length(λ)
        println("λ$i = $(λ[i])")
        println("v$i = $(v[:,i])")
    end
    return(nothing)
end

if length(ARGS) < 1
    @time lya()
else
    @time test_run()
end
