# Script to construct figure 3
using Assembly
using Plots
using JLD
using DataFrames
using StatsPlots
using StatsBase
using Statistics
using LaTeXStrings
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

function figure3(Rls::Array{Int64,1},Rus::Array{Int64,1},syns::Array{Bool,1},ens::Array{String,1},Ni::Int64,Nr::Int64)
    # Check if all these vectors are the same length
    if length(Rls) != length(Rus) || length(Rls) != length(syns) || length(Rls) != length(ens)
        error("length of vectors doesn't match")
    end
    println("Compiled!")
    # Count number of parameter sets
    Ns = length(Rls)
    # Container to store number of survivors
    svs = zeros(Int64,Ns,Nr)
    # Container to store metabolite diversity
    mbs = zeros(Int64,Ns,Nr)
    # Container to store dissipation
    dsp = zeros(Float64,Ns,Nr)
    # Preallocate labels
    lbs = Array{String}(undef,Ns)
    # Loop over parameter sets
    for i = 1:Ns
        for j = 1:Nr
            # Read in relevant files
            pfile = "Data/$(Rls[i])-$(Rus[i])$(syns[i])$(Ni)$(ens[i])/RedParasReacs$(Rls[i])-$(Rus[i])Syn$(syns[i])Run$(j)Ns$(Ni).jld"
            if ~isfile(pfile)
                error("parameter set $(i) run $(j) is missing a parameter file")
            end
            ofile = "Data/$(Rls[i])-$(Rus[i])$(syns[i])$(Ni)$(ens[i])/RedOutputReacs$(Rls[i])-$(Rus[i])Syn$(syns[i])Run$(j)Ns$(Ni).jld"
            if ~isfile(ofile)
                error("parameter set $(i) run $(j) is missing an output file")
            end
            # Only want final parameter set
            ps = load(pfile,"ps")
            inf_out = load(ofile,"inf_out")
            # Save number of survivors
            svs[i,j] = ps.N
            # Loop over metabolites to find those with non-zero concentrations
            cm = 0 # Set up counter
            mm = 0 # Lowest metabolite
            for k = 2:ps.M
                if inf_out[ps.N+k] > 0.0
                    # Increment counter
                    cm += 1
                    # And save new minimum metabolite
                    mm = k
                end
            end
            # Save results to vector
            mbs[i,j] = cm
            # Find and save dissipation
            dsp[i,j] = dissipation(ps,ps.mics,inf_out)
        end
        # Make labels here
        if ens[i] == "h"
            if syns[i] == true
                lbs[i] = "high-true"
            else
                lbs[i] = "high-false"
            end
        elseif ens[i] == "i"
            if syns[i] == true
                lbs[i] = "inter-true"
            else
                lbs[i] = "inter-false"
            end
        else
            if syns[i] == true
                lbs[i] = "low-true"
            else
                lbs[i] = "low-false"
            end
        end
    end
    # Make empty containers to store data as 1D vector
    tl = Array{String,1}(undef,Ns*Nr)
    tsv = Array{Int64,1}(undef,Ns*Nr)
    tdv = Array{Int64,1}(undef,Ns*Nr)
    enl = Array{String,1}(undef,Ns*Nr)
    eps = Array{Float64,1}(undef,Ns*Nr)
    # Fill out with data
    for i = 1:Ns
        for j = 1:Nr
            tl[(i-1)*Nr+j] = lbs[i]
            tsv[(i-1)*Nr+j] = svs[i,j]
            tdv[(i-1)*Nr+j] = mbs[i,j]
            enl[(i-1)*Nr+j] = ens[i]
            eps[(i-1)*Nr+j] = dsp[i,j]
        end
    end
    # Collect everything into one data frame
    survivors = DataFrame(PSet=tl,ns=tsv,sdv=tdv,en=enl,ent=eps)
    # Find corresponding indices of reordered labels
    pos = zeros(Float64,Ns)
    # Find posistions for mean and std
    for i = 1:Ns
        p = 0.5
        if syns[i] == true
            p += 1.0
        end
        if ens[i] == "i"
            p += 2.4
        elseif ens[i] == "l"
            p += 4.8
        end
        pos[i] = p
    end
    # Make latex label
    JKs = L"JK^{-1}s^{-1}"
    # Setup plotting
    pyplot()
    theme(:wong2,dpi=200)
    wongc = get_color_palette(wong_palette,57)
    # Want to do the plotting here
    p1 = plot(title="Ecosystem diversity",ylabel="Number of surviving strains")
    @df survivors violin!(p1,:PSet,:ns,linewidth=0,label="",color=wongc[2],group=:en)
    @df survivors boxplot!(p1,:PSet,:ns,color=wongc[4],fillalpha=0.75,linewidth=2,label="",group=:en)
    # Plot means
    for i = 1:Ns
        # Calculate mean
        mn = mean(svs[i,:])
        sdn = std(svs[i,:])
        scatter!(p1,[pos[i]],[mn],yerror=[sdn],label="",shape=:star5,color=wongc[5],ms=10,msc=wongc[5])
    end
    # Calculate minimum and maximum diversities
    maxd = convert(Float64,maximum(svs))
    mind = convert(Float64,maximum(svs))
    # Add annotation
    px, py = annpos([-1.0],[maxd; mind])
    annotate!(px,py,text("B",17,:black))
    savefig(p1,"Output/Fig3/Diversity.png")
    p2 = plot(title="Substrate diversification",ylabel="Number of substrates")
    @df survivors violin!(p2,:PSet,:sdv,linewidth=0,label="",color=wongc[2],group=:en)
    @df survivors boxplot!(p2,:PSet,:sdv,color=wongc[4],fillalpha=0.75,linewidth=2,label="",group=:en)
    # Plot means
    for i = 1:Ns
        # Calculate mean
        mn = mean(mbs[i,:])
        sdn = std(mbs[i,:])
        scatter!(p2,[pos[i]],[mn],yerror=[sdn],label="",shape=:star5,color=wongc[5],ms=10,msc=wongc[5])
    end
    # Calculate minimum and maximum substrate diversities
    maxs = convert(Float64,maximum(mbs))
    mins = convert(Float64,maximum(mbs))
    # Add annotation
    px, py = annpos([-1.0],[maxs; mins])
    annotate!(px,py,text("C",17,:black))
    savefig(p2,"Output/Fig3/SubDiv.png")
    p3 = plot(title="Entropy production",ylabel="Ecosystem entropy production rate ($(JKs))",yaxis=:log10)
    @df survivors violin!(p3,:PSet,:ent,linewidth=0,label="",color=wongc[2],group=:en)
    @df survivors boxplot!(p3,:PSet,:ent,color=wongc[4],fillalpha=0.75,linewidth=2,label="",group=:en)
    # Plot means
    for i = 1:Ns
        # Calculate mean
        mn = mean(dsp[i,:])
        sdn = std(dsp[i,:])
        scatter!(p3,[pos[i]],[mn],yerror=[sdn],label="",shape=:star5,color=wongc[5],ms=10,msc=wongc[5])
    end
    # Calculate minimum and maximum entropies
    maxe = convert(Float64,maximum(dsp))
    mine = convert(Float64,maximum(dsp))
    # Add annotation
    px, py = annpos([-1.0],[maxe; mine])
    annotate!(px,py,text("D",17,:black))
    savefig(p3,"Output/Fig3/EntropyProduction.png")
    # NOW WANT TO MAKE PLOT OF STRAIN DIVERSITY WITH TIME
    p4 = plot(title="Placeholder")
    # Now want to make a plot incorperating all four previous plots
    pt = plot(p4,p2,p1,p3,layout=4,size=(1200,800))
    savefig(pt,"Output/Fig3/figure3.png")
    return(nothing)
end

# Hard code parameters here
l = [1,1,1,1,1,1]
u = [5,5,5,5,5,5]
s = [true,true,true,false,false,false]
e = ["l","i","h","l","i","h"]

@time figure3(l,u,s,e,250,250)
