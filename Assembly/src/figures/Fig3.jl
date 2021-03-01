# Script to construct figure 3
using Assembly
using Plots
using JLD
using DataFrames
using StatsPlots
using StatsBase
using Statistics
using LaTeXStrings
using ColorSchemes
using Distributions
import PyPlot

# Function do Welch's t test
function my_t_test(data1::Array,data2::Array)
    # Find sample sizes
    N1 = length(data1)
    N2 = length(data2)
    # Preallocate means and standard deviation
    sam_m = zeros(2)
    sam_sd = zeros(2)
    # Calculate means and standard deviations
    sam_m[1] = mean(data1)
    sam_sd[1] = std(data1)
    sam_m[2] = mean(data2)
    sam_sd[2] = std(data2)
    # Calculate t value
    t = (sam_m[1]-sam_m[2])/sqrt((sam_sd[1]^2)/N1+(sam_sd[2]^2)/N2)
    # Calculate ν value
    ν = (((sam_sd[1]^2)/N1+(sam_sd[2]^2)/N2)^2)/((sam_sd[1]^4)/((N1^2)*(N1-1))+(sam_sd[2]^4)/((N2^2)*(N2-1)))
    # make corresponding t distribution
    T = TDist(ν)
    # Two tailed test
    P = 2*(1-cdf(T,abs(t)))
    return(P)
end

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

# Function to make bar chart of the diversity loss
function divloss(Rl::Int64,Ru::Int64,syn::Bool,en::String,Ni::Int64,Nr::Int64)
    # Read in relevant files
    pfile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/ParasReacs$(Rl)-$(Ru)Syn$(syn)Run$(Nr)Ns$(Ni).jld"
    if ~isfile(pfile)
        error("run $(Nr) is missing a parameter file")
    end
    ofile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/OutputReacs$(Rl)-$(Ru)Syn$(syn)Run$(Nr)Ns$(Ni).jld"
    if ~isfile(ofile)
        error("run $(Nr) is missing an output file")
    end
    efile = "Data/$(Rl)-$(Ru)$(syn)$(Ni)$(en)/ExtinctReacs$(Rl)-$(Ru)Syn$(syn)Run$(Nr)Ns$(Ni).jld"
    if ~isfile(efile)
        error("run $(Nr) is missing an extinct file")
    end
    # Load everything that's potentially useful
    ps = load(pfile,"ps")
    C = load(ofile,"C")
    T = load(ofile,"T")
    out = load(ofile,"out")
    ded = load(efile,"ded")
    # Remake inital list of microbes
    ms = Array{MicrobeP,1}(undef,Ni)
    # Setup counter
    cnt = 0
    # Loop over all microbes
    for i = 1:Ni
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
    r_cat = zeros(Int64,Ni)
    # Loop over every microbe in here and catagorise
    for i = 1:Ni
        # Find index of main reaction
        _, ind = findmax(ms[i].ϕP)
        # Find reaction
        r = ps.reacs[ms[i].Reacs[ind]]
        # And store index
        r_cat[i] = r.Rct
    end
    # Choose number of time points to do
    Tp = 8
    # Find and save maximum time
    Tmax = T[end]
    # Make vector of the time points
    # IS THIS A SENSIBLE TIME RANGE???
    Ts = [collect(range(0.0,stop=Tmax/25.0,length=Tp-1)); Tmax]
    # Container to time varying diversity and abundance data
    abT = zeros(Float64,Tp,ps.M)
    # Loop over reactions
    for i = 1:ps.M
        # Find all indices for this reaction
        inds = findall(x->x==i,r_cat)
        for j = 1:Tp
            # Find first index past time point
            Tind = findfirst(x->x>=Ts[j],T)
            # Sum abundances of survivors
            if Tind != 1
                # Interpolate abundances
                dT = ((Ts[j]-T[Tind-1])*sum(C[Tind,inds]) + (T[Tind]-Ts[j])*sum(C[Tind-1,inds]))/(T[Tind]-T[Tind-1])
                # Then save
                abT[j,i] = dT
            else
                abT[j,i] = sum(C[Tind,inds])
            end
        end
    end
    # Rescale bars to 1
    for i = 1:Tp
        abT[i,:] = abT[i,:]/(sum(abT[i,:]))
    end
    # Create xlabels
    xs = fill("",Tp)
    for i = 1:Tp
        if Ts[i] != 0.0
            # Round time
            Tr = round(Ts[i],sigdigits=2)
            # Find power of ten
            p10 = floor(Int64,log10(Tr))
            # reduce T by this factor of 10
            rT = Tr/(10.0^p10)
            # Put together as label
            xs[i] = "$(rT)e$(p10)"
        else
            xs[i] = "0.0"
        end
    end
    # Load in colorschemes
    a = ColorSchemes.tab20.colors
    b = ColorSchemes.tab20b.colors
    # Empty array to store colors
    cls = Array{RGB{Float64},1}(undef,40)
    for i = 1:20
        cls[i] = a[i]
    end
    for i = 21:40
        cls[i] = b[i-20]
    end
    # Make my own color scheme
    sch = ColorScheme(cls)
    # Now setup plotting
    pyplot(dpi=200)
    p = groupedbar(abT,bar_position=:stack,label="",palette=sch)
    plot!(p,xticks=(1:Tp,xs),ylabel="Relative abundance",xlabel="Time (s)")
    plot!(p,title="Diversity with time")
    # Add annotation
    px, py = annpos([0.25; convert(Float64,Tp)],[1.0;0.0],0.10,0.05)
    annotate!(px,py,text("A",17,:black))
    savefig(p,"Output/Fig3/abT.png")
    return(p)
end

function figure3(Rls::Array{Int64,1},Rus::Array{Int64,1},syns::Array{Bool,1},ens::Array{String,1},
                Ni::Int64,Nr::Int64,dRl::Int64,dRu::Int64,dsyn::Bool,den::String,runN::Int64)
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
    spm = Array{Float64,1}(undef,Ns*Nr)
    # Fill out with data
    for i = 1:Ns
        for j = 1:Nr
            tl[(i-1)*Nr+j] = lbs[i]
            tsv[(i-1)*Nr+j] = svs[i,j]
            tdv[(i-1)*Nr+j] = mbs[i,j]
            enl[(i-1)*Nr+j] = ens[i]
            eps[(i-1)*Nr+j] = dsp[i,j]
            spm[(i-1)*Nr+j] = svs[i,j]/mbs[i,j]
        end
    end
    # Collect everything into one data frame
    survivors = DataFrame(PSet=tl,ns=tsv,sdv=tdv,en=enl,ent=eps,rat=spm)
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
    # Container to store mean + sd for each case
    msd = zeros(Ns)
    # Find variables to loop over
    vars = unique(ens)
    # Add a t test here
    for i = 1:length(vars)
        # find energy condition to look for
        en = vars[i]
        # Find indices to do comparisons over
        inds = findall(x->x==en,ens)
        # Calculate all three
        Pd = my_t_test(svs[inds[1],:],svs[inds[2],:])
        Ps = my_t_test(mbs[inds[1],:],mbs[inds[2],:])
        Pe = my_t_test(dsp[inds[1],:],dsp[inds[2],:])
        Pr = my_t_test(svs[inds[1],:]./mbs[inds[1],:],svs[inds[2],:]./mbs[inds[2],:])
        println("Looking at $(en) case")
        println("Diversity P value = $(Pd)")
        println("Substrate diversification P value = $(Ps)")
        println("Entropy production P value = $(Pe)")
        println("Ratio P value = $(Pr)")
    end
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
        # sdn = sem(svs[i,:])
        msd[i] = mn + sdn
        scatter!(p1,[pos[i]],[mn],yerror=[sdn],label="",shape=:star5,color=wongc[5],ms=10,msc=wongc[5])
    end
    # Calculate minimum and maximum diversities
    maxd = convert(Float64,maximum(svs))
    mind = convert(Float64,maximum(svs))
    # Add annotation
    px, py = annpos([0.2;7.0],[maxd; mind; msd],0.10,0.05)
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
        # sdn = sem(mbs[i,:])
        msd[i] = mn + sdn
        scatter!(p2,[pos[i]],[mn],yerror=[sdn],label="",shape=:star5,color=wongc[5],ms=10,msc=wongc[5])
    end
    # Calculate minimum and maximum substrate diversities
    maxs = convert(Float64,maximum(mbs))
    mins = convert(Float64,maximum(mbs))
    # Add annotation
    px, py = annpos([0.2;7.0],[maxs; mins; msd],0.10,0.05)
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
        # sdn = sem(dsp[i,:])
        msd[i] = mn + sdn
        scatter!(p3,[pos[i]],[mn],yerror=[sdn],label="",shape=:star5,color=wongc[5],ms=10,msc=wongc[5])
    end
    # Calculate minimum and maximum entropies
    maxe = convert(Float64,maximum(dsp))
    mine = convert(Float64,maximum(dsp))
    # Add annotation
    px, py = annpos([0.2;7.0],[maxe; mine; msd],0.10,0.2)
    annotate!(px,py,text("D",17,:black))
    savefig(p3,"Output/Fig3/EntropyProduction.png")
    # Run div loss function to make extra plot
    p4 = divloss(dRl,dRu,dsyn,den,Ni,runN)
    # Now want to make a plot incorperating all four previous plots
    pt = plot(p4,p2,p1,p3,layout=4,size=(1200,800))
    savefig(pt,"Output/Fig3/figure3.png")
    # Want to do the plotting here
    plot(title="Ratio",ylabel="Survivors per substrate")
    @df survivors violin!(:PSet,:rat,linewidth=0,label="",color=wongc[2],group=:en)
    @df survivors boxplot!(:PSet,:rat,color=wongc[4],fillalpha=0.75,linewidth=2,label="",group=:en)
    # Plot means
    for i = 1:Ns
        # Calculate mean
        mr = mean(svs[i,:]./mbs[i,:])
        sdr = std(svs[i,:]./mbs[i,:])
        # sdr = sem(svs[i,:]./mbs[i,:])
        scatter!([pos[i]],[mr],yerror=[sdr],label="",shape=:star5,color=wongc[5],ms=10,msc=wongc[5])
    end
    savefig("Output/Fig3/Ratio.png")
    return(nothing)
end

# DELETE THIS WHEN DONE
# @time divloss(1,5,true,"i",250,89)

# Hard code parameters here
l = [1,1,1,1,1,1]
u = [5,5,5,5,5,5]
s = [true,true,true,false,false,false]
e = ["l","i","h","l","i","h"]

@time figure3(l,u,s,e,250,250,1,5,true,"i",89)
