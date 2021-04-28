# Script to construct figure 4
using Assembly
using Plots
using JLD
using StatsBase
using Plots.PlotMeasures
using LaTeXStrings
import PyPlot

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

function figure4(Rl::Int64,Ru::Int64,syns::Array{Bool,1},ens::Array{String,1},rps::Int64,
                Ni::Int64,fT::Float64,sT::Float64)
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
    # Read in first data file (chose l)
    ofile = "Data/$(Rl)-$(Ru)$(syns[1])$(Ni)l/OutputReacs$(Rl)-$(Ru)Syn$(syns[1])Run1Ns$(Ni).jld"
    if ~isfile(ofile)
        error("run 1 energy supply l is missing an output file")
    end
    pfile = "Data/$(Rl)-$(Ru)$(syns[1])$(Ni)l/ParasReacs$(Rl)-$(Ru)Syn$(syns[1])Run1Ns$(Ni).jld"
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
            if isfile("Output/Fig4/$(syns[l])$(ens[j]).jld")
                println("Reading $(ens[j])-$(syns[l]) data")
                # If it does then load the data in
                tfile = "Output/Fig4/$(syns[l])$(ens[j]).jld"
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
                    pfile = "Data/$(Rl)-$(Ru)$(syns[l])$(Ni)$(ens[j])/ParasReacs$(Rl)-$(Ru)Syn$(syns[l])Run$(i)Ns$(Ni).jld"
                    if ~isfile(pfile)
                        error("parameter set $(i) run $(j) is missing a parameter file")
                    end
                    ofile = "Data/$(Rl)-$(Ru)$(syns[l])$(Ni)$(ens[j])/OutputReacs$(Rl)-$(Ru)Syn$(syns[l])Run$(i)Ns$(Ni).jld"
                    if ~isfile(ofile)
                        error("run $(i) energy supply $(ens[j]) is missing an output file")
                    end
                    efile = "Data/$(Rl)-$(Ru)$(syns[l])$(Ni)$(ens[j])/ExtinctReacs$(Rl)-$(Ru)Syn$(syns[l])Run$(i)Ns$(Ni).jld"
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
                # Write out data here
                jldopen("Output/Fig4/$(syns[l])$(ens[j]).jld","w") do file
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
    savefig(p1,"Output/Fig4/DvTime$(Rl)-$(Ru).png")
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
    savefig(p2,"Output/Fig4/SubDvTime$(Rl)-$(Ru).png")
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
    savefig(p3,"Output/Fig4/Efficency.png")
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
    savefig(p4,"Output/Fig4/GrowthRate.png")
    # Plot all graphs as a single figure
    pt = plot(p1,p4,p2,p3,layout=4,size=(1200,800),margin=7.5mm)
    savefig(pt,"Output/Fig4/figure4.png")
    return(nothing)
end

@time figure4(1,5,[true,false],["l","h"],250,250,1e-2,0.0)
