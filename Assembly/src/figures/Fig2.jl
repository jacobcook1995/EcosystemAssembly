# Script to construct figure 2
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

function figure2(Rl::Int64,Ru::Int64,syn::Bool,Nr::Int64,Ns::Int64,en::String,Tf::Float64,rps::Int64)
    println("Compiled!")
    # Set intial number of concentrations, and preallocate the final ones
    nmi = ones(rps)
    nmf = zeros(rps)
    # Loop over all repeats to find substrate diversification
    for i = 1:rps
        # Read in relevant files
        pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/RedParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ns).jld"
        if ~isfile(pfile)
            error("run $(i) is missing a parameter file")
        end
        ofile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/RedOutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(i)Ns$(Ns).jld"
        if ~isfile(ofile)
            error("run $(Nr) is missing an output file")
        end
        # Load required data
        ps = load(pfile,"ps")
        inf_out = load(ofile,"inf_out")
        nmf[i] = count(x->x>0.0,inf_out[(ps.N+1):(ps.N+ps.M)])
    end
    # Read in specific files needed for the dynamics
    pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/ParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(Nr)Ns$(Ns).jld"
    if ~isfile(pfile)
        error("run $(Nr) is missing a parameter file")
    end
    ofile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/OutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(Nr)Ns$(Ns).jld"
    if ~isfile(ofile)
        error("run $(Nr) is missing an output file")
    end
    efile = "Data/$(Rl)-$(Ru)$(syn)$(Ns)$(en)/ExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(Nr)Ns$(Ns).jld"
    if ~isfile(efile)
        error("run $(Nr) is missing an extinct file")
    end
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
    # Strains to display
    ds = 25
    # Check if this is enough to show all survivors
    if ds < ps.N
        error("number of survivors higher than number of strains shown")
    end
    # Find indicies of strains I want to show
    is = zeros(Int64,ds)
    for i = 1:ps.N
        is[i] = findfirst(x->x==out[i],C[end,:])
    end
    # Fill in remaining non-survivors
    for i = (ps.N+1):ds
        # Find strains not already included
        dff = setdiff(1:Ns,is[1:i-1])
        is[i] = dff[1]
    end
    # Set suitable threshold for viable usage of metabolite 3
    tsh3 = 1e-4
    # Find and save time where metabolite three has first crossed threshold
    ind3 = findfirst(x->x>=tsh3,C[:,Ns+3])
    T3 = T[ind3]
    # Now move onto plotting
    pyplot()
    theme(:wong2,dpi=200)
    # Find maximum number of substrates (convert to integer)
    mS = convert(Int64,maximum(nmf))
    # make appropriate bins
    rbins = range(-0.25,stop=mS+0.75,length=mS+2)
    # make histograms
    ph = plot(xlabel="Number of substrates")
    histogram!(ph,nmi,color=:black,label="",bins=rbins)
    histogram!(ph,nmf,color=:red,label="",bins=rbins)
    savefig(ph,"Output/Fig2/hist.png")
    # Plot all the populations
    p1 = plot(title="Population dynamics",yaxis=:log10,xlabel="Time (s)",ylabel="Population (# cells)")
    # Store max and min C values
    maxC = zeros(length(is))
    minC = zeros(length(is))
    c = 0
    for i = is
        c += 1
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0) .& (T .<= Tend)
        plot!(p1,T[inds],C[inds,i],label="")
        # Store max and min C values in range
        maxC[c] = maximum(C[inds,i])
        minC[c] = minimum(C[inds,i])
    end
    # Add annotation
    px, py = annpos([0.0; Tend],[maxC; minC])
    # Different because log10 scale used
    annotate!(px,py,text("A",17,:black),δx=0.10,δy=0.125)
    vline!(p1,[T3],color=:red,style=:dash,label="")
    savefig(p1,"Output/Fig2/pops.png")
    # Now plot concentrations
    p2 = plot(title="Substrate diversification",xlabel="Time (s)",ylabel="Concentration (moles)")
    # Store max and min C values
    maxC = zeros(length(is))
    minC = zeros(length(is))
    c = 0
    for i = (Ns+1):(Ns+ps.M)
        c += 1
        # Find and eliminate points after end time
        inds = (T .<= Tend)
        # Can't switch theme but can switch pallete to avoid repeated colors
        if (i-Ns) >= 20
            plot!(p2[1],T[inds],C[inds,i],label="",palette=:darktest)
        else
            plot!(p2[1],T[inds],C[inds,i],label="")
        end
        # Store max and min C values in range
        maxC[c] = maximum(C[inds,i])
        minC[c] = minimum(C[inds,i])
    end
    # Add annotation
    px, py = annpos([0.0; Tend],[maxC; minC])
    annotate!(px,py,text("B",17,:black))
    vline!(p2[1],[T3],color=:red,style=:dash,label="")
    # Define box for inset here
    box = (1,bbox(0.7,0.25,0.275,0.275,:bottom,:left))
    histogram!(p2,nmi,color=:black,label="Initial",bins=rbins,inset_subplots=box,subplot=2)
    histogram!(p2[2],nmf,color=:red,label="Final",bins=rbins,xlabel="Number of substrates")
    plot!(p2[2],guidefontsize=6,legendfontsize=6,tickfontsize=4)
    savefig(p2,"Output/Fig2/concs.png")
    # Now plot proteome fraction
    p3 = plot(title="Proteome dynamics",xlabel="Time (s)",ylabel="Ribosome fraction")
    # Store max and min C values
    maxC = zeros(length(is))
    minC = zeros(length(is))
    c = 0
    for i = is.+2*Ns.+ps.M
        c += 1
        # Find and eliminate points after end time, remove points where strain is dead
        inds = (C[:,i-2*Ns-ps.M] .> 0) .& (T .<= Tend)
        plot!(p3,T[inds],C[inds,i],label="")
        # Store max and min C values in range
        maxC[c] = maximum(C[inds,i])
        minC[c] = minimum(C[inds,i])
    end
    # Add annotation
    px, py = annpos([0.0; Tend],[maxC; minC])
    annotate!(px,py,text("D",17,:black))
    vline!(p3,[T3],color=:red,style=:dash,label="")
    savefig(p3,"Output/Fig2/fracs.png")
    # container to store entropy production
    ep = zeros(length(T))
    # Calculate entropy production at each step
    for i = 1:length(T)
        # Calculate entropy production at each step
        ep[i] = dissipation(ps,ms,C[i,:])
    end
    p4 = plot(title="Entropy production",xlabel="Time (s)",ylabel="Entropy production (J/K per s)")
    # Find and eliminate points after end time
    inds = (T .<= Tend)
    plot!(p4,T[inds],ep[inds],label="")
    # Add annotation
    px, py = annpos([0.0; Tend],ep[inds])
    annotate!(px,py,text("C",17,:black))
    vline!(p4,[T3],color=:red,style=:dash,label="")
    savefig(p4,"Output/Fig2/entp.png")
    # Now want to make a plot incorperating all four previous plots
    pt = plot(p1,p4,p2,p3,layout=4,size=(1200,800))
    savefig(pt,"Output/Fig2/figure2.png")
    return(nothing)
end

@time figure2(1,5,true,61,250,"i",0.01,250)
