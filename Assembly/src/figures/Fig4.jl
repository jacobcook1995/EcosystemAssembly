# Script to construct figure 4
using Assembly
using Plots
using JLD
using StatsBase
using Plots.PlotMeasures
import PyPlot

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
    dv = zeros(ips,rps,L)
    tab = zeros(ips,rps,L)
    # Set threshold for substrate being properly diversified
    tsh = 1e-7
    # Calculate length
    le = length(ens)
    # Loop over repeats
    for l = 1:length(syns)
        for j = 1:length(ens)
            for i = 1:rps
                # Read in relevant files
                ofile = "Data/$(Rl)-$(Ru)$(syns[l])$(Ni)$(ens[j])/OutputReacs$(Rl)-$(Ru)Syn$(syns[l])Run$(i)Ns$(Ni).jld"
                if ~isfile(ofile)
                    error("run $(i) energy supply $(ens[j]) is missing an output file")
                end
                # Load full dynamics
                C = load(ofile,"C")
                T = load(ofile,"T")
                # Loop over the time points
                for k = 1:ips
                    # Find first time point greater than or equal to one were looking for
                    ind = findfirst(x->x>=Ts[k],T)
                    # If time points are equal just save number of survivors
                    if T[ind] == Ts[k]
                        svs[k,i,(j-1)*le+l] = count(x->x>0.0,C[ind,1:Ni])
                        dv[k,i,(j-1)*le+l] = count(x->x>0.0,C[ind,(Ni+1):(Ni+ps.M-1)])
                        tab[k,i,(j-1)*le+l] = sum(C[ind,1:Ni])
                    else
                        # Otherwise need to (linearly) interpolate
                        dT = (T[ind]-Ts[k])/(T[ind]-T[ind-1])
                        svs[k,i,(j-1)*le+l] = (1-dT)*count(x->x>0.0,C[ind,1:Ni]) + dT*count(x->x>0.0,C[ind-1,1:Ni])
                        dv[k,i,(j-1)*le+l] = (1-dT)*count(x->x>0.0,C[ind,(Ni+1):(Ni+ps.M-1)]) + dT*count(x->x>0.0,C[ind-1,(Ni+1):(Ni+ps.M-1)])
                        tab[k,i,(j-1)*le+l] = (1-dT)*sum(C[ind,1:Ni]) + dT*sum(C[ind,1:Ni])
                    end
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
    # Find mean and standard errors of survivor numbers and substrate diversification
    for i = 1:ips
        for j = 1:length(ens)
            for k = 1:length(syns)
                msvs[i,(j-1)*le+k] = mean(svs[i,:,(j-1)*le+k])
                sdsvs[i,(j-1)*le+k] = sem(svs[i,:,(j-1)*le+k])
                mdv[i,(j-1)*le+k] = mean(dv[i,:,(j-1)*le+k])
                sddv[i,(j-1)*le+k] = sem(dv[i,:,(j-1)*le+k])
                mta[i,(j-1)*le+k] = mean(tab[i,:,(j-1)*le+k])
                sdta[i,(j-1)*le+k] = sem(tab[i,:,(j-1)*le+k])
            end
        end
    end
    # Set up plotting
    pyplot()
    # Set a color-blind friendly palette
    theme(:wong2,dpi=200)
    # NEED TO FIND AND PLOT SOME COMMEN LINE HERE!!!!
    # Make labels
    lb = ["low-true" "low-false" "high-true" "high-false"]
    # Plot survivors
    p1 = plot(title="Diversity with time",xlabel="Time",ylabel="Number of surviving strains")
    plot!(p1,Ts,msvs[:,1:L],ribbon=sdsvs[:,1:L],labels=lb,legend=:right)
    # Twin the xaxis and then plot substrate diversity
    plot!(p1,twinx(),Ts,mdv[:,1:L],ribbon=sddv[:,1:L],labels="",ylabel="Number of substrates")
    # Increase left margin
    plot!(p1,right_margin=20.0mm)
    savefig(p1,"Output/Fig4/DvTime$(Rl)-$(Ru).png")
    # Plot total abundances
    p2 = plot(title="Ecosystem abundance with time",xlabel="Time",ylabel="Total abundance")
    plot!(p2,yaxis=:log10)
    for i = 1:L
        plot!(p2,Ts,mta[:,i],ribbon=sdta[:,i],label=lb[i])
    end
    savefig(p2,"Output/Fig4/TotalAbTime$(Rl)-$(Ru).png")
    return(nothing)
end

@time figure4(1,5,[true,false],["l","h"],250,250,2.5e-2,0.0)
