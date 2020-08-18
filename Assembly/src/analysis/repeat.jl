# A script to do repeated runs of the Inhibition model, adding new microbes
# in if they can successfully invade
using Assembly
using JLD
using SymPy

# function to run simulation to steady state
function stead(ps::InhibParameters,Tmax::Float64,pop::Array{Float64,1},conc::Array{Float64,1})
    # Preallocate vector of forces
    F = Array{Sym,1}(undef,ps.N+ps.M)
    # Find force vector using function
    F = Force(ps,F)
    # Preallocate output
    C = Array{Float64,2}(undef,0,ps.N+ps.M)
    T = Array{Float64,1}(undef,0)
    # Now start doing runs
    std = false
    t = Tmax
    c = 0
    while std == false
        c += 1
        # Run the simulation
        Ci, Ti = inhib_simulate(ps,t,pop,conc)
        # Use final simulation results to find local forces
        f = nForce(F,Ci[end,:],ps)
        # Check maximum eigenvalues
        if maximum(abs.(f)) <= 1.0e-5 || c == 100
            # Either end loop
            println("Steady state found after $c steps")
            if c == 100
                println("MASSIVE PROBLEM HERE!")
            end
            std = true
        else
            # Otherwise update initial condiditions
            pop = Ci[end,1:ps.N]
            conc = Ci[end,ps.N+1:end]
            # Switch to using small time step to smooth out numerical oscillations
            if c == 5
                t = 5.0
            elseif c == 10
                t = 1.0
            elseif c == 15
                t = 0.1
            end
            if c != 0 && c % 10 == 0
                println("$c steps taken without finding steady state")
            end
        end
        # Regardless of outcome add simulation to output
        if length(T) > 0
            T = cat(T,Ti.+T[end],dims=1)
        else
            T = cat(T,Ti,dims=1)
        end
        C = cat(C,Ci[:,:],dims=1)
    end
    return(C,T)
end

# function to check if microbe can actually grow
function check(ps::InhibParameters,pop::Array{Float64,1},conc::Array{Float64,1})
    # Preallocate Jacobian
    F = Array{Sym,1}(undef,ps.N+ps.M)
    # Find Jacobian using function
    F = Force(ps,F)
    # Use initial conditions to find local forces
    f = nForce(F,[pop;1.0;conc],ps)
    # Check if any microbe grows
    if maximum(f[ps.N]) > 1.0e-5
        return(true)
    else
        return(false)
    end
end

# function to repeatedly add microbes
function repeat()
    println("Successfully compiled.")
    if length(ARGS) == 0
        error("NEED TO PROVIDE NAME FOR OUTPUT DATA.")
    end
    # Make new parameter set if one isn't provided
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
    # Now attempt to make the parameter set
    vld = false
    ps = 0 # This put it in the right scope
    while vld == false
        ps = initialise(N,M,O,mR,sdR,mq,sdq,mK,sdK,mk,sdk)
        # Only accept parameter set if at least one reaction leaves metabolite 1
        if any((ps.reacs.↦:Rct) .== 1)
            vld = true
        end
    end
    # NEED A BETTER NAMING CUSTOM HERE
    # save this parameter set
    jldopen("Paras/ps$(ARGS[1]).jld","w") do file
        write(file,"ps",ps)
    end
    # Then simulate to steady state
    Tmax = 100.0
    # Initialise vectors of concentrations and populations
    pop = ones(ps.N)
    conc = zeros(ps.M) # No chemical to begin with
    # Now find first steady state
    C1, T1 = stead(ps,Tmax,pop,conc)
    # Check if microbe survives
    if C1[end,1] > 0.0
        surv = true
        exT = 0
    else
        surv = false
        exT = 1
    end
    # Then use relevant information to construct metacommunity
    MD = make_MetaCom([make_MicData(ps.mics[1],true,0,surv,exT)])
    # Make vector of addition times for microbes
    t = [0.0]
    # Want to try to add 10 microbes
    for i = 1:250
        pop = C1[end,1:ps.N]
        conc = C1[end,ps.N+1:end]
        # Use function to construct random microbe based on parameter set
        mic = initialise_mic(ps,mR,sdR,mq,sdq,mK,sdK,mk,sdk)
        # Add microbe to temporary parameter set
        ps2 = add_Microbe(ps,mic)
        # Now use function to check if this microbe will grow
        chk = check(ps2,pop,conc)
        if chk == true
            t = cat(t,T1[end],dims=1)
            println("Microbe $i should grow.")
            # Update parameter set
            ps = ps2
            Cn, Tn = stead(ps,Tmax,[pop;1.0],conc)
            # Can just cat T straight on
            T1 = cat(T1,Tn.+T1[end],dims=1)
            # Need to update C1 so that it can be cat'ed to
            CT = zeros(size(C1,1),size(C1,2)+1)
            CT[:,1:ps.N-1] = C1[:,1:ps.N-1]
            CT[:,ps.N] .= 0.0
            CT[:,ps.N+1:end] = C1[:,ps.N:end]
            # Then finally cat CT here
            C1 = cat(CT,Cn,dims=1)
            # Check if microbe survives
            if C1[end,length(t)] > 0.0
                md = make_MicData(mic,true,i,true,0)
            else
                println("Microbe doesn't survive after all")
                md = make_MicData(mic,true,i,false,length(t))
            end
            # Add this successful invader to the data
            MD = add_MetaCom(MD,md)
            # Finally check if any other microbes have been gone extinct
            for j = 1:length(t)
                if MD.data[j].surv == true
                    if C1[end,j] == 0.0
                        # Use function to add an extinction to the data set
                        MD = extinct_MetaCom(MD,j,length(t))
                    end
                end
            end
        else
            println("Microbe $i could not grow.")
            # Add this failed invader to the data
            md = make_MicData(mic,false,i,false,length(t))
            MD = add_MetaCom(MD,md)
        end
    end
    # save this metacommunity data
    jldopen("Output/MD$(ARGS[1]).jld","w") do file
        write(file,"MD",MD)
    end
    # WHAT DATA SHOULD BE OUTPUTTED AND SAVED?
    jldopen("Output/Dyn$(ARGS[1]).jld","w") do file
        write(file,"C",C1)
        write(file,"T",T1)
        write(file,"t",t)
    end
    return(nothing)
end

# function repeat_test()
#     # Load in previous data
#     ps = load("Temp/Paras/ps$(ARGS[1]).jld","ps")
#     MD = load("Temp/Paras/MD$(ARGS[1]).jld","MD")
#     # Then simulate to steady state
#     Tmax = 100.0
#     # Initialise vectors of concentrations and populations
#     pop = ones(ps.N)
#     conc = zeros(ps.M) # No chemical to begin with
#     # Now find first steady state
#     C1, T1 = stead(ps,Tmax,pop,conc)
#     # Make vector of addition times for microbes
#     t = [0.0]
#     # Run for first five microbes
#     for i = 1:10
#         pop = C1[end,1:ps.N]
#         conc = C1[end,ps.N+1:end]
#         # Use function to construct random microbe based on parameter set
#         mic = MD.data[i+1].mic
#         # Add microbe to temporary parameter set
#         ps2 = add_Microbe(ps,mic)
#         # Now use function to check if this microbe will grow
#         chk = check(ps2,pop,conc)
#         if chk == true
#             t = cat(t,T1[end],dims=1)
#             println("Microbe $i should grow.")
#             # Update parameter set
#             ps = ps2
#             Cn, Tn = stead(ps,Tmax,[pop;1.0],conc)
#             # Can just cat T straight on
#             T1 = cat(T1,Tn.+T1[end],dims=1)
#             # Need to update C1 so that it can be cat'ed to
#             CT = zeros(size(C1,1),size(C1,2)+1)
#             CT[:,1:ps.N-1] = C1[:,1:ps.N-1]
#             CT[:,ps.N] .= 0.0
#             CT[:,ps.N+1:end] = C1[:,ps.N:end]
#             # Then finally cat CT here
#             C1 = cat(CT,Cn,dims=1)
#         else
#             println("Microbe $i could not grow.")
#         end
#     end
#     return(nothing)
# end

# Plotting function
function plot_repeat()
    if length(ARGS) == 0
        error("NEED TO PROVIDE NAME FOR INPUT DATA.")
    end
    # First read in parameter set
    ps = load("Paras/ps$(ARGS[1]).jld","ps")
    # Make reaction matrix
    rM = zeros(Int64,ps.M,ps.M)
    for i = 1:length(ps.reacs)
        rM[ps.reacs[i].Prd,ps.reacs[i].Rct] += 1
    end
    pyplot(dpi=200)
    heatmap(1:size(rM,1),1:size(rM,2),rM,colorbar=:none)
    plot!(xlabel="Reactant",ylabel="Product")
    plot!(0:100,0:100,color=:white,label="")
    savefig("Output/RandReacs$(ARGS[1]).png")
    # Now read in dynamics data
    C = load("Output/Dyn$(ARGS[1]).jld","C")
    T = load("Output/Dyn$(ARGS[1]).jld","T")
    t = load("Output/Dyn$(ARGS[1]).jld","t")
    # Need to recalculate N
    N = ps.N + length(t) - 1
    # Now plot the dynamics
    plot(T,C[:,1:N],label="",xlabel="Time",ylabel="Population")
    vline!(t,color=:red,label="",style=:dot)
    savefig("Output/Pops$(ARGS[1]).png")
    plot(T,C[:,N+1:end],label="")
    vline!(t,color=:red,label="",style=:dot)
    savefig("Output/Concs$(ARGS[1]).png")
    return(nothing)
end

@time repeat()
