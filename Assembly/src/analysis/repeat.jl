# A script to do repeated runs of the Inhibition model, adding new microbes
# in if they can successfully invade
using Assembly
using Plots
import PyPlot
using JLD
using SymPy

# function to run simulation to steady state
function stead(ps,Tmax,pop,conc)
    # Preallocate Jacobian
    J = Array{Sym,2}(undef,ps.N+ps.M,ps.N+ps.M)
    # Find Jacobian using function
    J = Jacobian(ps,J)
    # Preallocate output
    C = Array{Float64,2}(undef,ps.N+ps.M,0)
    T = Array{Float64,1}(undef,0)
    # Now start doing runs
    std = false
    c = 1
    while std == false
        Ci, Ti = inhib_simulate(ps,Tmax,pop,conc)
        # Use final simulation results to find local Lyapunov exponents
        λ, _ = Lyapunov(J,Ci[end,:],ps)
        # Check maximum eigenvalues
        if maximum(real(λ)) < 0.0
            # Either end loop
            std = true
        else
            # Otherwise update initial condiditions
            pop = Ci[end,1:ps.N]
            conc = Ci[end:ps.N+1:end]
            println("Failed")
        end
        # Regardless of outcome add simulation to output
        if c == 1
            T = Ti
            C = Ci
        else
            T = cat(T,Ti,dims=1)
            C = cat(C,Ci,dims=2)
        end
    end
    return(C,T)
end

# function to check if microbe can actually grow
function check(ps::InhibParameters,pop::Array{Float64,1},conc::Array{Float64,1})
    # Preallocate Jacobian
    J = Array{Sym,2}(undef,ps.N+ps.M,ps.N+ps.M)
    # Find Jacobian using function
    J = Jacobian(ps,J)
    # Use initial conditions to find eigenvales
    λ, _ = Lyapunov(J,[pop;0.01;conc],ps)
    # Check largest eigenvalue
    if maximum(real(λ)) > 0.0
        return(true)
    else
        return(false)
    end
end

# function to repeatedly add microbes
function repeat()
    println("Successfully compiled.")
    # Initially make 1 microbe
    N = 1
    M = 100
    O = 2*M
    mR = 10.0
    sdR = 1.0
    mq = 1.0
    sdq = 0.1
    mK  = 0.1
    sdK = 0.01
    mk = 10.0
    sdk = 1.0
    # Now make the parameter set
    ps = initialise(N,M,O,mR,sdR,mq,sdq,mK,sdK,mk,sdk)
    # save this parameter set
    jldopen("Temp/Paras/psR.jld","w") do file
        write(file,"ps",ps)
    end
    # Then simulate to steady state
    Tmax = 100.0
    # Initialise vectors of concentrations and populations
    pop = ones(ps.N)
    conc = zeros(ps.M) # No chemical to begin with
    # Now find first steady state
    C1, T1 = stead(ps,Tmax,pop,conc)
    # Want to try to add 10 microbes
    for i = 1:10
        # Use function to construct random microbe based on parameter set
        mic = initialise_mic(ps,mR,sdR,mq,sdq,mK,sdK,mk,sdk)
        # Add microbe to temporary parameter set
        ps2 = add_Microbe(ps,mic)
        # Now use function to check if this microbe will grow
        chk = check(ps2,pop,conc)
        if chk == true
            println("Microbe $i grows")
            # Update parameter set
            ps = ps2
            Cn, Tn = stead(ps,Tmax,[pop;0.01],conc)
            # Can just cat T straight on
            T1 = cat(T1,Tn.+T1[end],dims=1)
            # Need to update C1 so that it can be cat'ed to
            CT = zeros(size(C1,dims=1),size(C1,dims=1)+1)
            CT[:,1:ps.N-1] = C1[:,1:ps.N-1]
            CT[:,ps.N] .= 0.0
            CT[:,ps.N+1:end] = C1[:,ps.N:end]
            # Then finally cat CT here
            C1 = cat(CT,Cn,dims=2)
        else
            println("Microbe $i could not grow")
        end
    end
    pyplot(dpi=200)
    plot(T1,C1[:,1:ps.N],label="")
    savefig("Output/TestPop.png")
    plot(T1,C1[:,ps.N+1:end],label="")
    savefig("Output/TestConc.png")
return(nothing)
end

@time repeat()
