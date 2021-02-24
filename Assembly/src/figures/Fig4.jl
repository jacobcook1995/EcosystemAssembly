# Script to construct figure 4
using Assembly
using Plots
using JLD
using StatsBase
using Plots.PlotMeasures
import PyPlot

# Function to calculate fraction of thermodynamically unfeasible reactions
function therm_unf(pop::Array{Float64,1},conc::Array{Float64,1},ms::Array{MicrobeP,1},ps::FullParameters)
    # Find indices of surving populations
    inds = findall(x->x>0.0,pop)
    # Setup two seperate counters
    RT = 0 # Total number of reactions
    Ruf = 0 # Number of thermodynamically unfeasible reactions
    # Loop over indices
    for i = inds
        # Find number of reactions to loop over
        R = ms[i].R
        # Add this to the total
        RT += R
        # Loop over reactions
        for j = 1:R
            # Find corresponding reaction
            r = ps.reacs[ms[i].Reacs[j]]
            # Check that substrate concentration is non-zero
            if conc[r.Rct] != 0.0
                # Find theta value
                θt = θ_smooth(conc[r.Rct],conc[r.Prd],ps.T,ms[i].η[j],r.ΔG0)
                # Check if reaction is unfeasible
                if θt > 0.99
                    Ruf += 1
                end
            end
        end
    end
    # Return as a proportion unfeasible
    return(Ruf/RT)
end


# Function to find proportion of reactions with no substrate
function no_sub(pop::Array{Float64,1},conc::Array{Float64,1},ms::Array{MicrobeP,1},ps::FullParameters)
    # Find indices of surving populations
    inds = findall(x->x>0.0,pop)
    # Setup two seperate counters
    RT = 0 # Total number of reactions
    Rns = 0 # Number of reactions with no substrate
    # Loop over indices
    for i = inds
        # Find number of reactions to loop over
        R = ms[i].R
        # Add this to the total
        RT += R
        # Loop over reactions
        for j = 1:R
            # Find corresponding reaction
            r = ps.reacs[ms[i].Reacs[j]]
            # Check if substrate concentration is zero
            if conc[r.Rct] == 0.0
                Rns += 1
            end
        end
    end
    # Return as a proportion unfeasible
    return(Rns/RT)
end

function av_eff(pop::Array{Float64,1},conc::Array{Float64,1},ms::Array{MicrobeP,1},ps::FullParameters)
    # Define mimimum product to substrate ratio (to calculate) the efficency
    mr = 1e-2
    # Find indices of surving populations
    inds = findall(x->x>0.0,pop)
    # Set initial efficency value
    efT = 0.0
    # Loop over indices
    for i = inds
        # Find number of reactions to loop over
        R = ms[i].R
        # Loop over reactions
        for j = 1:R
            # Find corresponding reaction
            r = ps.reacs[ms[i].Reacs[j]]
            # Calculate proportional effiency
            eff = -(ms[i].η[j]*ΔGATP)/(r.ΔG0+Rgas*ps.T*log(mr))
            # Add this to total weighted by population, and metabolic protein fraction
            efT += pop[i]*ms[i].ϕP[j]*eff
        end
    end
    # Average across populations
    avef = efT/sum(pop)
    return(avef)
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
    # New containers, for non-viable reactions and efficencies
    ufR = zeros(ips,rps,L)
    nsR = zeros(ips,rps,L)
    efs = zeros(ips,rps,L)
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
                ufR[:,:,(j-1)*le+l] = load(tfile,"ufR")
                nsR[:,:,(j-1)*le+l] = load(tfile,"nsR")
                efs[:,:,(j-1)*le+l] = load(tfile,"efs")
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
                            ufR[k,i,(j-1)*le+l] = therm_unf(C[ind,1:Ni],C[ind,(Ni+1):(Ni+ps.M)],ms,ps)
                            nsR[k,i,(j-1)*le+l] = no_sub(C[ind,1:Ni],C[ind,(Ni+1):(Ni+ps.M)],ms,ps)
                            efs[k,i,(j-1)*le+l] = av_eff(C[ind,1:Ni],C[ind,(Ni+1):(Ni+ps.M)],ms,ps)
                        else
                            # Otherwise need to (linearly) interpolate
                            dT = (T[ind]-Ts[k])/(T[ind]-T[ind-1])
                            svs[k,i,(j-1)*le+l] = (1-dT)*count(x->x>0.0,C[ind,1:Ni]) + dT*count(x->x>0.0,C[ind-1,1:Ni])
                            via[k,i,(j-1)*le+l] = (1-dT)*count(x->x>=popt,C[ind,1:Ni]) + dT*count(x->x>=popt,C[ind-1,1:Ni])
                            dv[k,i,(j-1)*le+l] = (1-dT)*count(x->x>0.0,C[ind,(Ni+1):(Ni+ps.M-1)]) + dT*count(x->x>0.0,C[ind-1,(Ni+1):(Ni+ps.M-1)])
                            tab[k,i,(j-1)*le+l] = (1-dT)*sum(C[ind,1:Ni]) + dT*sum(C[ind,1:Ni])
                            ufR[k,i,(j-1)*le+l] = (1-dT)*therm_unf(C[ind,1:Ni],C[ind,(Ni+1):(Ni+ps.M)],ms,ps)
                            ufR[k,i,(j-1)*le+l] += dT*therm_unf(C[ind-1,1:Ni],C[ind-1,(Ni+1):(Ni+ps.M)],ms,ps)
                            nsR[k,i,(j-1)*le+l] = (1-dT)*no_sub(C[ind,1:Ni],C[ind,(Ni+1):(Ni+ps.M)],ms,ps)
                            nsR[k,i,(j-1)*le+l] += dT*no_sub(C[ind-1,1:Ni],C[ind-1,(Ni+1):(Ni+ps.M)],ms,ps)
                            efs[k,i,(j-1)*le+l] = (1-dT)*av_eff(C[ind,1:Ni],C[ind,(Ni+1):(Ni+ps.M)],ms,ps)
                            efs[k,i,(j-1)*le+l] += dT*av_eff(C[ind-1,1:Ni],C[ind-1,(Ni+1):(Ni+ps.M)],ms,ps)
                        end
                    end
                end
                # Write out data here
                jldopen("Output/Fig4/$(syns[l])$(ens[j]).jld","w") do file
                    write(file,"svs",svs[:,:,(j-1)*le+l])
                    write(file,"via",via[:,:,(j-1)*le+l])
                    write(file,"dv",dv[:,:,(j-1)*le+l])
                    write(file,"tab",tab[:,:,(j-1)*le+l])
                    write(file,"ufR",ufR[:,:,(j-1)*le+l])
                    write(file,"nsR",nsR[:,:,(j-1)*le+l])
                    write(file,"efs",efs[:,:,(j-1)*le+l])
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
    mufR = zeros(ips,L)
    sdufR = zeros(ips,L)
    mnsR = zeros(ips,L)
    sdnsR = zeros(ips,L)
    mimR = zeros(ips,L)
    sdimR = zeros(ips,L)
    mefs = zeros(ips,L)
    sdefs = zeros(ips,L)
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
                mufR[i,(j-1)*le+k] = mean(ufR[i,:,(j-1)*le+k])
                sdufR[i,(j-1)*le+k] = sem(ufR[i,:,(j-1)*le+k])
                mnsR[i,(j-1)*le+k] = mean(nsR[i,:,(j-1)*le+k])
                sdnsR[i,(j-1)*le+k] = sem(nsR[i,:,(j-1)*le+k])
                mimR[i,(j-1)*le+k] = mean(nsR[i,:,(j-1)*le+k].+ufR[i,:,(j-1)*le+k])
                sdimR[i,(j-1)*le+k] = sem(nsR[i,:,(j-1)*le+k].+ufR[i,:,(j-1)*le+k])
                mefs[i,(j-1)*le+k] = mean(efs[i,:,(j-1)*le+k])
                sdefs[i,(j-1)*le+k] = sem(efs[i,:,(j-1)*le+k])
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
    # Set up plotting
    pyplot()
    # Set a color-blind friendly palette
    theme(:wong2,dpi=200)
    wongc = get_color_palette(wong_palette,57)
    # Make labels
    lb = ["low-true" "low-false" "high-true" "high-false"]
    # Plot survivors
    p1 = plot(title="Diversity with time",xlabel="Time",ylabel="Number of surviving strains")
    vline!(p1,vTs,linestyle=:dash,label="")
    plot!(p1,Ts,msvs[:,1:L],ribbon=sdsvs[:,1:L],labels=lb,legend=:right)
    # Twin the xaxis and then plot substrate diversity
    plot!(twinx(),Ts,mdv[:,1:L],ribbon=sddv[:,1:L],labels="",ylabel="Number of substrates",palette=wongc[2:5])
    # Increase left margin
    plot!(p1,right_margin=20.0mm)
    # Add annotation
    px, py = annpos(Ts,[5.0,250.0],0.15,0.05)
    annotate!(p1,px,py,text("A",17,:black))
    savefig(p1,"Output/Fig4/DvTime$(Rl)-$(Ru).png")
    # Plot total abundances
    p2 = plot(title="Ecosystem abundance with time",xlabel="Time",ylabel="Total abundance")
    plot!(p2,yaxis=:log10,legend=:right)
    vline!(p2,vTs,linestyle=:dash,label="")
    for i = 1:L
        plot!(p2,Ts,mta[:,i],ribbon=sdta[:,i],label=lb[i])
    end
    # Add annotation
    maxab = vec(mta.+sdta)
    px, py = annpos(Ts,maxab,0.15,0.2)
    annotate!(p2,px,py,text("B",17,:black))
    savefig(p2,"Output/Fig4/TotalAbTime$(Rl)-$(Ru).png")
    # PLACEHOLDER GRAPHS, NEED TO DECIDE WHICH OF THESE TO KEEP
    p3 = plot(title="Unfeasible reactions",xlabel="Time",ylabel="Proportion of unfeasible reactions")
    plot!(p3,Ts,mufR[:,1:L],ribbon=sdufR[:,1:L],labels=lb)
    savefig(p3,"Output/Fig4/Unfeasible.png")
    p4 = plot(title="Reactions with no substrate",xlabel="Time",ylabel="Proportion of reactions without substrate")
    plot!(p4,Ts,mnsR[:,1:L],ribbon=sdnsR[:,1:L],labels=lb)
    savefig(p4,"Output/Fig4/Nosub.png")
    p5 = plot(title="Impossible reactions",xlabel="Time",ylabel="Proportion of impossible reactions")
    plot!(p5,Ts,mimR[:,1:L],ribbon=sdimR[:,1:L],labels=lb)
    savefig(p5,"Output/Fig4/ImpossibleReactions.png")
    p6 = plot(title="Average efficency with time",xlabel="Time",ylabel="Efficency of reactions")
    plot!(p6,Ts,mefs[:,1:L],ribbon=sdefs[:,1:L],labels=lb)
    savefig(p6,"Output/Fig4/Efficency.png")
    # Plot all graphs as a single figure
    pt = plot(p1,p3,p2,p4,layout=(2,2),size=(1200,800),margin=15.0mm)
    savefig(pt,"Output/Fig4/figure4.png")
    return(nothing)
end

@time figure4(1,5,[true,false],["l","h"],250,250,2.5e-2,0.0)
