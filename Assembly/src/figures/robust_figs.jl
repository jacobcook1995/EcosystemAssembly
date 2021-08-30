# Script to make robustness figures
using Assembly
using Plots
using JLD
using StatsPlots
using StatsBase
using Statistics
using Plots.PlotMeasures
using LaTeXStrings
using ColorSchemes
using KernelDensity
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

# function to find average efficency
function av_eff(pop::Array{Float64,1},conc::Array{Float64,1},ms::Array{MicrobeP,1},ps::FullParameters)
    # Define mimimum product to substrate ratio (to calculate) the efficency
    mr = 1e-2
    # Find indices of surving populations
    inds = findall(x->x>0.0,pop)
    # Set initial efficency value
    efT = 0.0
    # Loop over indices
    for i = inds
        # Loop over number of reactions
        for j = 1:ms[i].R
            # Calculate proportional effiency
            eff = -(ms[i].η[j]*ΔGATP)/(ps.reacs[ms[i].Reacs[j]].ΔG0+Rgas*ps.T*log(mr))
            # Add this to total weighted by population, and metabolic protein fraction
            efT += pop[i]*ms[i].ϕP[j]*eff
        end
    end
    # Average across populations
    avef = efT/sum(pop)
    return(avef)
end

# funtion to find average growth rate
function av_λ(pop::Array{Float64,1},as::Array{Float64,1},ϕRs::Array{Float64,1},ms::Array{MicrobeP,1})
    # Find indices of surving populations
    inds = findall(x->x>0.0,pop)
    # weighted total growth rate (starts at zero)
    λT = 0.0
    # Loop over survivors
    for i = inds
        # Find growth rate for this particular strain
        λt = λs(as[i],ϕRs[i],ms[i])
        # Weight by population and add to total
        λT += pop[i]*λt
    end
    # Divide weighted growth rate by total abundance
    λT /= sum(pop)
    return(λT)
end

function rob_figure3(Rls::Array{Int64,1},Rus::Array{Int64,1},syns::Array{Bool,1},ens::Array{String,1},
                Ni::Int64,Nr::Int64,opt::Int64)
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
    # Container to store functional diversities
    fdv = zeros(Int64,Ns,Nr)
    # Find title of option
    title_op = options_titles(opt)
    # Loop over parameter sets
    for i = 1:Ns
        for j = 1:Nr
            # Read in relevant files
            pfile = "Paras/$(title_op)/$(Rls[i])-$(Rus[i])$(syns[i])$(ens[i])/RedParasReacs$(Rls[i])-$(Rus[i])Syn$(syns[i])Run$(j)Ns$(Ni).jld"
            if ~isfile(pfile)
                error("parameter set $(i) run $(j) is missing a parameter file")
            end
            ofile = "Data/$(title_op)/$(Rls[i])-$(Rus[i])$(syns[i])$(ens[i])/RedOutputReacs$(Rls[i])-$(Rus[i])Syn$(syns[i])Run$(j)Ns$(Ni).jld"
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
            # Preallocate vector to store if function is present
            fnc = fill(false,ps.M)
            # Loop over every surviving microbe
            for k = 1:ps.N
                # Find index of main reaction
                _, ind = findmax(ps.mics[k].ϕP)
                # Find corresponding reactant number
                Rn = ps.reacs[ps.mics[k].Reacs[ind]].Rct
                # Update vector to show that function Rn is present
                if fnc[Rn] == false
                    fnc[Rn] = true
                end
            end
            # Count number of functions
            fdv[i,j] = count(fnc)
        end
    end
    # Check if directory exists and if not make it
    if ~isdir("Output/$(title_op)")
        mkdir("Output/$(title_op)")
    end
    # Check if sub-directory exists and if not make it
    if ~isdir("Output/$(title_op)/Fig3")
        mkdir("Output/$(title_op)/Fig3")
    end
    # Setup plotting
    pyplot()
    theme(:wong2,dpi=300)
    wongc = wong2_palette()
    # Find corresponding indices of reordered labels
    pos = zeros(Float64,Ns)
    c = Array{RGBA,1}(undef,Ns)
    # Find posistions for mean and std
    for i = 1:Ns
        p = 1.0
        if syns[i] == true
            p += 0.5
            # Set color
            c[i] = wongc[1]
        else
            # Choose nice grey as "light" color
            c[i] = RGB(([170, 170, 170] / 255)...)
        end
        if ens[i] == "l"
            p += 1.5
        end
        pos[i] = p
    end
    # Make latex label
    JKs = L"JK^{-1}s^{-1}"
    # Container to store mean + sd for each case
    msd = zeros(Ns)
    # Want to do the plotting here
    p1 = plot(ylabel="Number of surviving strains",xlim=(0.5,3.5),xlabel="Energy supply")
    plot!(p1,xticks=([1.25,2.75],["high","low"]))
    # Plot means
    for i = 1:Ns
        # Calculate mean
        mn = mean(svs[i,:])
        # Calculate 99% confidence interval
        sdn = sem(svs[i,:])*2.576
        scatter!(p1,[pos[i]],[mn],yerror=[sdn],label="",color=c[i],ms=6,msc=c[i])
    end
    # Add bracket for significance plot
    plot!(p1,[2.5,3.0],[5.0,5.0],color=:black,label="")
    plot!(p1,[2.5,2.5],[4.6,5.01],color=:black,label="")
    plot!(p1,[3.0,3.0],[4.6,5.01],color=:black,label="")
    # Then add star above the bracket
    scatter!(p1,[2.75],[5.25],color=:black,shape=:star6,label="")
    savefig(p1,"Output/$(title_op)/Fig3/Diversity.png")
    p3 = plot(ylabel="Entropy production rate ($(JKs))",yaxis=:log10)
    plot!(p3,xlim=(0.5,3.5),xticks=([1.25,2.75],["high","low"]),xlabel="Energy supply")
    # Plot means
    for i = 1:Ns
        # Calculate mean
        mn = mean(dsp[i,:])
        # Calculate 99% confidence interval
        sdn = sem(dsp[i,:])*2.576
        scatter!(p3,[pos[i]],[mn],yerror=[sdn],label="",color=c[i],ms=6,msc=c[i])
    end
    # Add bracket for significance plot
    plot!(p3,[2.5,3.0],[6.25e-2,6.25e-2],color=:black,label="")
    plot!(p3,[2.5,2.5],[5.4e-2,6.26e-2],color=:black,label="")
    plot!(p3,[3.0,3.0],[5.4e-2,6.26e-2],color=:black,label="")
    # Then add star above the bracket
    scatter!(p3,[2.75],[6.9e-2],color=:black,shape=:star6,label="")
    savefig(p3,"Output/$(title_op)/Fig3/EntropyProduction.png")
    # Want to do the plotting here
    p2 = plot(ylabel="Survivors per substrate",xlim=(0.5,3.5))
    plot!(p2,xticks=([1.25,2.75],["high","low"]),xlabel="Energy supply")
    # Plot means
    for i = 1:Ns
        # Calculate mean
        mn = mean(svs[i,:]./mbs[i,:])
        # Calculate 99% confidence interval
        sdn = sem(svs[i,:]./mbs[i,:])*2.576
        scatter!(p2,[pos[i]],[mn],yerror=[sdn],label="",color=c[i],ms=6,msc=c[i])
    end
    savefig(p2,"Output/$(title_op)/Fig3/Ratio.png")
    # Want to do the plotting here
    p4 = plot(ylabel="Number of surviving functional groups",xlim=(0.5,3.5))
    plot!(p4,xticks=([1.25,2.75],["high","low"]),xlabel="Energy supply",legend=:right)
    # Plot means
    for i = 1:Ns
        # Calculate mean
        mn = mean(fdv[i,:])
        # Calculate 99% confidence interval
        sdn = sem(fdv[i,:])*2.576
        # Save for annotation
        msd[i] = mn + sdn
        # Label empty
        lb = ""
        # Unless energy is high
        if ens[i] == "h"
            if syns[i] == true
                lb = "Reversible"
            else
                lb = "M–M"
            end
        end
        scatter!(p4,[pos[i]],[mn],yerror=[sdn],label=lb,color=c[i],ms=6,msc=c[i])
    end
    # Add annotation
    px, py = annpos([0.5;3.5],msd,0.5,-0.01)
    annotate!(px,py,text("B",17,:black))
    # Add bracket for significance plot
    plot!(p4,[2.5,3.0],[4.0,4.0],color=:black,label="")
    plot!(p4,[2.5,2.5],[3.7,4.01],color=:black,label="")
    plot!(p4,[3.0,3.0],[3.7,4.01],color=:black,label="")
    # Then add star above the bracket
    scatter!(p4,[2.75],[4.2],color=:black,shape=:star6,label="")
    savefig(p4,"Output/$(title_op)/Fig3/FuncDiv.png")
    # Combine all three plots into a single one
    pc = plot(p4,p1,p2,p3,layout=(1,4),size=(800,400),guidefontsize=13,legendfontsize=8,tickfontsize=11)
    savefig(pc,"Output/$(title_op)/Fig3/figure3.png")
    return(nothing)
end

function rob_figure4(Rl::Int64,Ru::Int64,syns::Array{Bool,1},ens::Array{String,1},rps::Int64,
                Ni::Int64,fT::Float64,sT::Float64,opt::Int64)
    println("Compiled!")
    # Check I haven't provided stupid times
    if fT <= 0.0 || fT > 1.0
        error("final time can't be greater than 1, or less than zero")
    end
    if sT >= fT || sT < 0.0
        error("starting time can't be equal to final time, or less than zero")
    end
    # Choosing to sample a thousand points for now
    ips = 1000
    # Find title of the option
    title_op = options_titles(opt)
    # Read in first data file (chose l)
    ofile = "Data/$(title_op)/$(Rl)-$(Ru)$(syns[1])l/OutputReacs$(Rl)-$(Ru)Syn$(syns[1])Run1Ns$(Ni).jld"
    if ~isfile(ofile)
        error("run 1 energy supply l is missing an output file")
    end
    pfile = "Paras/$(title_op)/$(Rl)-$(Ru)$(syns[1])l/ParasReacs$(Rl)-$(Ru)Syn$(syns[1])Run1Ns$(Ni).jld"
    if ~isfile(ofile)
        error("run 1 energy supply l is missing an parameter file")
    end
    ps = load(pfile,"ps")
    T = load(ofile,"T")
    # And extract maximum time
    Tmax = T[end]
    # Make vector of times to check at
    Ts = collect(range(Tmax*sT,Tmax*fT,length=ips))
    # Preallocate number of survivors, substrate diversities, total abundance
    L = length(ens)*length(syns)
    svs = zeros(ips,rps,L)
    via = zeros(ips,rps,L)
    dv = zeros(ips,rps,L)
    tab = zeros(ips,rps,L)
    # New containers, for efficencies and growth rates
    efs = zeros(ips,rps,L)
    λts = zeros(ips,rps,L)
    # Set threshold for substrate being properly diversified
    tsh = 1e-7
    # Threshold meaningful/viable population
    popt = 1e-5
    # Calculate length
    le = length(ens)
    # Loop over repeats
    for l = 1:length(syns)
        for j = 1:length(ens)
            # Check if data has already been generated
            if isfile("Output/$(title_op)/Fig4/$(syns[l])$(ens[j]).jld")
                println("Reading $(ens[j])-$(syns[l]) data")
                # If it does then load the data in
                tfile = "Output/$(title_op)/Fig4/$(syns[l])$(ens[j]).jld"
                svs[:,:,(j-1)*le+l] = load(tfile,"svs")
                via[:,:,(j-1)*le+l] = load(tfile,"via")
                dv[:,:,(j-1)*le+l] = load(tfile,"dv")
                tab[:,:,(j-1)*le+l] = load(tfile,"tab")
                efs[:,:,(j-1)*le+l] = load(tfile,"efs")
                λts[:,:,(j-1)*le+l] = load(tfile,"λts")
            else
                println("Generating $(ens[j])-$(syns[l]) data")
                for i = 1:rps
                    # Read in relevant files
                    pfile = "Paras/$(title_op)/$(Rl)-$(Ru)$(syns[l])$(ens[j])/ParasReacs$(Rl)-$(Ru)Syn$(syns[l])Run$(i)Ns$(Ni).jld"
                    if ~isfile(pfile)
                        error("parameter set $(i) run $(j) is missing a parameter file")
                    end
                    ofile = "Data/$(title_op)/$(Rl)-$(Ru)$(syns[l])$(ens[j])/OutputReacs$(Rl)-$(Ru)Syn$(syns[l])Run$(i)Ns$(Ni).jld"
                    if ~isfile(ofile)
                        error("run $(i) energy supply $(ens[j]) is missing an output file")
                    end
                    efile = "Data/$(title_op)/$(Rl)-$(Ru)$(syns[l])$(ens[j])/ExtinctReacs$(Rl)-$(Ru)Syn$(syns[l])Run$(i)Ns$(Ni).jld"
                    if ~isfile(efile)
                        error("run $(Nr) is missing an extinct file")
                    end
                    # Load full dynamics
                    C = load(ofile,"C")
                    T = load(ofile,"T")
                    out = load(ofile,"out")
                    ps = load(pfile,"ps")
                    ded = load(efile,"ded")
                    # Make new vector of microbes
                    ms = Array{MicrobeP,1}(undef,Ni)
                    # Setup counter
                    cnt = 0
                    # Loop over all microbes
                    for k = 1:Ni
                        # Check if it is a survivor
                        if C[end,k] != 0.0 && C[end,k] ∈ out
                            # If it is find and save it
                            ind = findfirst(x->x==C[end,k],out)
                            ms[k] = ps.mics[ind]
                        else
                            # Update counter
                            cnt += 1
                            # Use next element from ded vector
                            ms[k] = ded[cnt]
                        end
                    end
                    # Loop over the time points
                    for k = 1:ips
                        # Find first time point greater than or equal to one were looking for
                        ind = findfirst(x->x>=Ts[k],T)
                        # If time points are equal just save number of survivors
                        if T[ind] == Ts[k]
                            svs[k,i,(j-1)*le+l] = count(x->x>0.0,C[ind,1:Ni])
                            via[k,i,(j-1)*le+l] = count(x->x>=popt,C[ind,1:Ni])
                            dv[k,i,(j-1)*le+l] = count(x->x>0.0,C[ind,(Ni+1):(Ni+ps.M-1)])
                            tab[k,i,(j-1)*le+l] = sum(C[ind,1:Ni])
                            efs[k,i,(j-1)*le+l] = av_eff(C[ind,1:Ni],C[ind,(Ni+1):(Ni+ps.M)],ms,ps)
                            λts[k,i,(j-1)*le+l] = av_λ(C[ind,1:Ni],C[ind,(Ni+ps.M+1):(2*Ni+ps.M)],C[ind,(2*Ni+ps.M+1):end],ms)
                        else
                            # Otherwise need to (linearly) interpolate
                            dT = (T[ind]-Ts[k])/(T[ind]-T[ind-1])
                            svs[k,i,(j-1)*le+l] = (1-dT)*count(x->x>0.0,C[ind,1:Ni]) + dT*count(x->x>0.0,C[ind-1,1:Ni])
                            via[k,i,(j-1)*le+l] = (1-dT)*count(x->x>=popt,C[ind,1:Ni]) + dT*count(x->x>=popt,C[ind-1,1:Ni])
                            dv[k,i,(j-1)*le+l] = (1-dT)*count(x->x>0.0,C[ind,(Ni+1):(Ni+ps.M-1)]) + dT*count(x->x>0.0,C[ind-1,(Ni+1):(Ni+ps.M-1)])
                            tab[k,i,(j-1)*le+l] = (1-dT)*sum(C[ind,1:Ni]) + dT*sum(C[ind,1:Ni])
                            efs[k,i,(j-1)*le+l] = (1-dT)*av_eff(C[ind,1:Ni],C[ind,(Ni+1):(Ni+ps.M)],ms,ps)
                            efs[k,i,(j-1)*le+l] += dT*av_eff(C[ind-1,1:Ni],C[ind-1,(Ni+1):(Ni+ps.M)],ms,ps)
                            λts[k,i,(j-1)*le+l] = (1-dT)*av_λ(C[ind,1:Ni],C[ind,(Ni+ps.M+1):(2*Ni+ps.M)],C[ind,(2*Ni+ps.M+1):end],ms)
                            λts[k,i,(j-1)*le+l] += dT*av_λ(C[ind-1,1:Ni],C[ind-1,(Ni+ps.M+1):(2*Ni+ps.M)],C[ind-1,(2*Ni+ps.M+1):end],ms)
                        end
                    end
                end
                # Check if directory exists and if not make it
                if ~isdir("Output/$(title_op)")
                    mkdir("Output/$(title_op)")
                end
                # Check if sub-directory exists and if not make it
                if ~isdir("Output/$(title_op)/Fig4")
                    mkdir("Output/$(title_op)/Fig4")
                end
                # Write out data here
                jldopen("Output/$(title_op)/Fig4/$(syns[l])$(ens[j]).jld","w") do file
                    write(file,"svs",svs[:,:,(j-1)*le+l])
                    write(file,"via",via[:,:,(j-1)*le+l])
                    write(file,"dv",dv[:,:,(j-1)*le+l])
                    write(file,"tab",tab[:,:,(j-1)*le+l])
                    write(file,"efs",efs[:,:,(j-1)*le+l])
                    write(file,"λts",λts[:,:,(j-1)*le+l])
                end
            end
        end
    end
    # Preallocate means and sds
    msvs = zeros(ips,L)
    sdsvs = zeros(ips,L)
    mdv = zeros(ips,L)
    sddv = zeros(ips,L)
    mta = zeros(ips,L)
    sdta = zeros(ips,L)
    mvia = zeros(ips,L)
    sdvia = zeros(ips,L)
    mefs = zeros(ips,L)
    sdefs = zeros(ips,L)
    mλts = zeros(ips,L)
    sdλts = zeros(ips,L)
    # Find mean and standard errors of survivor numbers and substrate diversification
    for i = 1:ips
        for j = 1:length(ens)
            for k = 1:length(syns)
                msvs[i,(j-1)*le+k] = mean(svs[i,:,(j-1)*le+k])
                sdsvs[i,(j-1)*le+k] = sem(svs[i,:,(j-1)*le+k])
                mvia[i,(j-1)*le+k] = mean(via[i,:,(j-1)*le+k])
                sdvia[i,(j-1)*le+k] = sem(via[i,:,(j-1)*le+k])
                mdv[i,(j-1)*le+k] = mean(dv[i,:,(j-1)*le+k])
                sddv[i,(j-1)*le+k] = sem(dv[i,:,(j-1)*le+k])
                mta[i,(j-1)*le+k] = mean(tab[i,:,(j-1)*le+k])
                sdta[i,(j-1)*le+k] = sem(tab[i,:,(j-1)*le+k])
                mefs[i,(j-1)*le+k] = mean(efs[i,:,(j-1)*le+k])
                sdefs[i,(j-1)*le+k] = sem(efs[i,:,(j-1)*le+k])
                mλts[i,(j-1)*le+k] = mean(λts[i,:,(j-1)*le+k])
                sdλts[i,(j-1)*le+k] = sem(λts[i,:,(j-1)*le+k])
            end
        end
    end
    # Preallocate viability times
    vTs = zeros(L)
    # Calculate viability times
    for i = 1:length(ens)
        for j = 1:length(syns)
            # find index of first time with 75% below threshold
            indT = findfirst(x->x<0.75*Ni,mvia[:,(i-1)*le+j])
            # Then save corresponding time
            vTs[(i-1)*le+j] = Ts[indT]
        end
    end
    # Set line width (plots defaults to 1)
    wdt = 1.5
    # Set up plotting
    pyplot(dpi=300,legendfontsize=10,tickfontsize=10,guidefontsize=12)
    # Making my own pallette
    pl = Array{RGBA,1}(undef,L)
    # Define the 4 colors I want to show
    pl[1] = RGB(([83,51,237] / 255)...) # dark-blue
    pl[2] = RGB(([137,196,244] / 255)...) # light-blue
    pl[3] = RGB(([207,0,15] / 255)...) # dark-red
    pl[4] = RGB(([241,169,160] / 255)...) # light-red
    # Make labels
    lb = ["low energy + reversible" "low energy + M-M" "high energy + reversible" "high energy + M-M"]
    # Plot survivors
    p1 = plot(title="Strain diversity",ylabel="Number of surviving strains")
    vline!(p1,vTs/1e6,linestyle=:dash,color=:black,label="")
    for i = 1:L
        plot!(p1,Ts/1e6,msvs[:,i],ribbon=sdsvs[:,i],label=lb[i],legend=:right,color=pl[i],lw=wdt)
    end
    # Add annotation
    px, py = annpos(Ts/1e6,[5.0,250.0],0.1,0.05)
    annotate!(p1,px,py,text("A",17,:black))
    savefig(p1,"Output/$(title_op)/Fig4/DvTime$(Rl)-$(Ru).png")
    # Make label
    s6 = L"10^6\;s"
    # plot substrate diversification
    p2 = plot(title="Substrate diversification",ylabel="Number of substrates",xlabel="Time ($(s6))",legend=false)
    vline!(p2,vTs/1e6,linestyle=:dash,color=:black,label="")
    for i = 1:L
        plot!(p2,Ts/1e6,mdv[:,i],ribbon=sddv[:,i],label=lb[i],color=pl[i],lw=wdt)
    end
    # Add annotation
    px, py = annpos(Ts/1e6,[0.0,23.0],0.1,0.05)
    annotate!(p2,px,py,text("B",17,:black))
    savefig(p2,"Output/$(title_op)/Fig4/SubDvTime$(Rl)-$(Ru).png")
    # Reduce end time for second two plots
    Tend = 1e6
    # Plot graph of efficencies
    p3 = plot(title="Average efficency",xlabel="Time ($(s6))",ylabel="Efficency of reactions",legend=false)
    for i = 1:L
        plot!(p3,Ts/1e6,mefs[:,i],ribbon=sdefs[:,i],label=lb[i],color=pl[i],lw=wdt)
    end
    # Add vertical line
    vline!(p3,vTs/1e6,linestyle=:dash,color=:black,label="")
    # Add annotation
    maxefs = vec(mefs.+sdefs)
    px, py = annpos([0.0,Tend/1e6],maxefs,0.1,0.1)
    annotate!(p3,px,py,text("D",17,:black))
    savefig(p3,"Output/$(title_op)/Fig4/Efficency.png")
    # Make label
    s1 = L"s^{-1}"
    # Plot graph of growth rates
    p4 = plot(title="Growth rate",ylabel="Growth rate ($(s1))",legend=false)
    for i = 1:L
        plot!(p4,Ts/1e6,mλts[:,i],ribbon=sdλts[:,i],label=lb[i],color=pl[i],lw=wdt)
    end
    # Add vertical line
    vline!(p4,vTs/1e6,linestyle=:dash,color=:black,label="")
    # Add annotation
    maxλts = vec(mλts.+sdλts)
    px, py = annpos([0.0,Tend/1e6],maxλts,0.1,0.05)
    annotate!(p4,px,py,text("C",17,:black))
    savefig(p4,"Output/$(title_op)/Fig4/GrowthRate.png")
    # Plot all graphs as a single figure
    pt = plot(p1,p4,p2,p3,layout=4,size=(1200,800),margin=7.5mm)
    savefig(pt,"Output/$(title_op)/Fig4/figure4.png")
    return(nothing)
end

# Hard code parameters here
l = [1,1,1,1]
u = [5,5,5,5]
s = [true,true,false,false]
e = ["l","h","l","h"]

@time rob_figure3(l,u,s,e,250,50)

@time rob_figure4(1,5,[true,false],["l","h"],50,250,1e-2,0.0,1)
