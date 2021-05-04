# NEED TO WORK OUT WHERE TO MOVE ALL THIS EVENTUALLY
using Assembly
using Plots
using JLD
import PyPlot

# function to calculate the dissipation for an assembled ecosystem
function dissipation(ps::FullParameters,ms::Array{MicrobeP,1},out::Array{Float64,1})
    # Set all elements of out less than zero to zero
    out[out.<0.0] .= 0.0
    # Define number of strains
    N = length(ms)
    # check that parameter set is sensible given the output
    if length(out) != ps.M + 3*N
        error("parameter set doesn't match output")
    end
    # Set dissipation to zero
    dsp = 0
    # Loop over number of strains
    for i = 1:N
        # Isolate this strain
        mic = ms[i]
        # Loop over reactions of this strain
        for j = 1:mic.R
            # Find appropriate reaction
            r = ps.reacs[mic.Reacs[j]]
            # If there's no product Gibbs free energy becomes infinite. Justified to ignore
            # this as if product hasn't built up reaction can't be happening to a sigificant degree
            if out[N+r.Prd] != 0.0
                # Find amount of energy that this reaction dissipates
                Fd = -(r.ΔG0 + Rgas*ps.T*log(out[N+r.Prd]/out[N+r.Rct]) + mic.η[j]*ΔGATP)
                # Find amount of enzyme E
                E = Eα(out[2*N+ps.M+i],mic,j)
                # Then find the rate that this reaction proceeds at
                q = qs(out[N+r.Rct],out[N+r.Prd],E,j,mic,ps.T,r)
                # Check if reaction actually occurs
                if q != 0.0
                    dsp += q*Fd*out[i]
                end
            end
        end
    end
    # Convert from molecule units to moles
    dsp /= NA
    return(dsp)
end

# Function to count the number of spikes in the entropy production trace
function ent_sp_cnt(Rl::Int64,Ru::Int64,syn::Bool,Ns::Int64,en::String,Tf::Float64,rps::Int64)
    println("Compiled!")
    # Loop over all repeats to find substrate diversification
    for i = 1:rps
        # Read in specific files needed for the dynamics
        pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/ParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ns).jld"
        if ~isfile(pfile)
            error("run $(Nr) is missing a parameter file")
        end
        ofile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/OutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ns).jld"
        if ~isfile(ofile)
            error("run $(Nr) is missing an output file")
        end
        efile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/ExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ns).jld"
        if ~isfile(efile)
            error("run $(Nr) is missing an extinct file")
        end
        # Read in relevant data
        ps = load(pfile,"ps")
        C = load(ofile,"C")
        T = load(ofile,"T")
        out = load(ofile,"out")
        ded = load(efile,"ded")
        # Make new vector of microbes
        ms = Array{MicrobeP,1}(undef,Ns)
        # Setup counter
        cnt = 0
        # Loop over all microbes
        for i = 1:Ns
            # Check if it is a survivor
            if C[end,i] != 0.0 && C[end,i] ∈ out
                # If it is find and save it
                ind = findfirst(x->x==C[end,i],out)
                ms[i] = ps.mics[ind]
            else
                # Update counter
                cnt += 1
                # Use next element from ded vector
                ms[i] = ded[cnt]
            end
        end
        # Find final time
        Tmax = T[end]
        # Find time to plot too
        Tend = Tmax*Tf
        # container to store entropy production
        ep = zeros(length(T))
        # Calculate entropy production at each step
        for i = 1:length(T)
            # Calculate entropy production at each step
            ep[i] = dissipation(ps,ms,C[i,:])
        end
        # container to store changes in entropy productions
        dep = zeros(length(T))
        # Find change between one point and the previous one
        for i = 1:length(T)
            if i > 1
                dep[i] = ep[i] - ep[i-1]
            else
                dep[i] = ep[i]
            end
        end
        # Find first peak, i.e. first time entropy production decays
        pind = findfirst(x->x<0.0,dep) - 1
        # Store in a vector
        pinds = [pind]
        # Counter for peaks
        pc = 1
        # Has the end of the trajectory been reached
        nd = false
        # Set inital offset
        off = pind
        # Now identify every peak
        while nd == false
            # Find first point where entropy production increases again
            try
                tind = findfirst(x->x>0.0,dep[(off+1):end]) - 1
                # Add trough to offset
                off += tind
            # Error arises because we're at the end of the vector
            catch e
                # Set offset to maximum
                off = length(T)
                # and finish loop
                nd = true
            end
            # Then find next peak based on this
            try
                pind = findfirst(x->x<0.0,dep[(off+1):end]) - 1
            # Error arises because we're at the end of the vector
            catch e
                # Set final peak
                pind = length(T) - off
                # and finish loop
                nd = true
            end
            # Check difference between peak and trough height
            df = (abs(ep[(pind+off)]) - abs(ep[(off)]))/(abs(ep[(pind+off)]))
            # If a greater than 1% change save
            if df > 0.01
                # Add new value to the vector
                pinds = cat(pinds,pind+off,dims=1)
                # Increment counter
                pc += 1
            end
            # Finally add 'peak' to the offset, regardless of if it's a true peak or not
            off += pind
        end
        println("Run $(i)")
        println(pc)
        # Now move onto plotting
        pyplot()
        theme(:wong2,dpi=300,guidefontsize=16,tickfontsize=14)
        # Set labels
        p4 = plot(xlabel="Time (s)",ylabel="Entropy production (J/K per s)")
        # Find and eliminate points after end time
        inds = (T .<= Tend)
        plot!(p4,T[inds],ep[inds],label="",ylim=(-0.01,Inf))
        for j = 1:length(pinds)
            if T[pinds[j]] <= Tend
                vline!(p4,[T[pinds[j]]],color=:red,linestyle=:dot,label="")
            end
        end
        savefig(p4,"Output/entp$(i).png")
    end
    return(nothing)
end

@time ent_sp_cnt(1,5,true,250,"i",1.0,250)
# @time ent_sp_cnt(1,5,true,250,"i",0.025,250)
