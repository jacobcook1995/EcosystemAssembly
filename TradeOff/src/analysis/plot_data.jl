# A script to read in and analyse model output data
using TradeOff
using JLD
using Plots
import PyPlot

# function to plot the trajectories
function plot_traj()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 3
        error("insufficent inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    rN = 0
    ims = 0
    sim_type = 0
    # Check that all arguments can be converted to integers
    try
        rN = parse(Int64,ARGS[1])
        ims = parse(Int64,ARGS[2])
        sim_type = parse(Int64,ARGS[3])
    catch e
            error("need to provide 3 integers")
    end
    println("Compiled")
    # Extract other simulation parameters from the function
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Read in appropriate files
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    ofile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/Run$(rN)Data$(ims)Ims.jld"
    if ~isfile(ofile)
        error("$(ims) immigrations run $(rN) is missing an output file")
    end
    # Read in relevant data
    ps = load(pfile,"ps")
    traj = load(ofile,"traj")
    T = load(ofile,"T")
    micd = load(ofile,"micd")
    its = load(ofile,"its")
    println("Data read in")
    # Find C from a function
    C = merge_data(ps,traj,T,micd,its)
    println("Data merged")
    # Check if directory exists and if not make it
    if ~isdir("Output/Plotsd=$(d)u=$(μrange)")
        mkdir("Output/Plotsd=$(d)u=$(μrange)")
    end
    # Find total number of strains
    totN = length(micd)
    # Find indices of extinct strains
    ext = isnan.(C[end,1:totN])
    # Invert to find survivors
    svr = .!ext
    pyplot(dpi=200)
    # Plot all the populations
    p1 = plot(yaxis=:log10,ylabel="Population (# cells)",ylims=(1e-5,Inf))
    for i = 1:totN
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,i] .> 0)
        plot!(p1,T[inds],C[inds,i],label="")
    end
    savefig(p1,"Output/Plotsd=$(d)u=$(μrange)/all_pops.png")
    # Plot all the concentrations
    p2 = plot(yaxis=:log10,ylabel="Concentration")#,ylims=(1e-15,Inf))
    for i = 1:ps.M
        # Find and eliminate zeros so that they can be plotted on a log plot
        inds = (C[:,totN+i] .> 0)
        plot!(p2,T[inds],C[inds,totN+i],label="")
    end
    savefig(p2,"Output/Plotsd=$(d)u=$(μrange)/all_concs.png")
    # Plot all the energy concentrations
    p3 = plot(ylabel="Energy Concentration")
    for i = 1:totN
        plot!(p3,T,C[:,totN+ps.M+i],label="")
    end
    savefig(p3,"Output/Plotsd=$(d)u=$(μrange)/all_as.png")
    # Plot all the ribosome fractions
    p4 = plot(ylabel="Ribosome fraction")
    for i = 1:totN
        plot!(p4,T,C[:,2*totN+ps.M+i],label="")
    end
    savefig(p4,"Output/Plotsd=$(d)u=$(μrange)/all_fracs.png")
    # Plot populations that survive to the end
    p1 = plot(yaxis=:log10,ylabel="Population (# cells)",ylims=(1e-5,Inf))
    for i = 1:totN
        if svr[i] == true
            # Find and eliminate zeros so that they can be plotted on a log plot
            inds = (C[:,i] .> 0)
            plot!(p1,T[inds],C[inds,i],label="")
        end
    end
    savefig(p1,"Output/Plotsd=$(d)u=$(μrange)/surv_pops.png")
    # Plot energy concentrations of populations that survive to the end
    p3 = plot(ylabel="Energy Concentration")
    for i = 1:totN
        if svr[i] == true
            plot!(p3,T,C[:,totN+ps.M+i],label="")
        end
    end
    savefig(p3,"Output/Plotsd=$(d)u=$(μrange)/surv_as.png")
    # Plot ribosome fractions of populations that survive to the end
    p4 = plot(ylabel="Ribosome fraction")
    for i = 1:totN
        if svr[i] == true
            plot!(p4,T,C[:,2*totN+ps.M+i],label="")
        end
    end
    savefig(p4,"Output/Plotsd=$(d)u=$(μrange)/surv_fracs.png")
    return(nothing)
end

# Function to load in and plot the averages
function plot_aves()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 2
        error("insufficent inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    ims = 0
    sim_type = 0
    # Check that all arguments can be converted to integers
    try
        ims = parse(Int64,ARGS[1])
        sim_type = parse(Int64,ARGS[2])
    catch e
        error("need to provide 2 integers")
    end
    println("Compiled")
    # Load in hardcoded simulation parameters
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Find file name to load in
    sfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)u=$(μrange)/RunStats$(ims)Ims.jld"
    # Check it actually exists
    if ~isfile(sfile)
        error("missing stats file for $(ims) immigrations simulations")
    end
    # Now load out the times, and number of trajectories
    times = load(sfile,"times")
    no_sims = load(sfile,"no_sims")
    no_via = load(sfile,"no_via")
    no_rs = load(sfile,"no_rs")
    # Load in averages
    mn_svt = load(sfile,"mn_svt")
    mn_tsvt = load(sfile,"mn_tsvt")
    mn_sbs = load(sfile,"mn_sbs")
    mn_Rs = load(sfile,"mn_Rs")
    mn_via_R = load(sfile,"mn_via_R")
    mn_ηs_R = load(sfile,"mn_ηs_R")
    mn_ωs_R = load(sfile,"mn_ωs_R")
    mn_kc_R = load(sfile,"mn_kc_R")
    mn_KS_R = load(sfile,"mn_KS_R")
    mn_kr_R = load(sfile,"mn_kr_R")
    mn_ηs = load(sfile,"mn_ηs")
    mn_via_η = load(sfile,"mn_via_η")
    mn_ωs = load(sfile,"mn_ωs")
    mn_via_ω = load(sfile,"mn_via_ω")
    mn_fr_ΔG = load(sfile,"mn_fr_ΔG")
    mn_via_a = load(sfile,"mn_via_a")
    mn_via_ϕR = load(sfile,"mn_via_ϕR")
    mn_kcs = load(sfile,"mn_kcs")
    mn_KSs = load(sfile,"mn_KSs")
    mn_krs = load(sfile,"mn_krs")
    mn_av_steps = load(sfile,"mn_av_steps")
    # Load in standard deviations
    sd_svt = load(sfile,"sd_svt")
    sd_tsvt = load(sfile,"sd_tsvt")
    sd_sbs = load(sfile,"sd_sbs")
    sd_Rs = load(sfile,"sd_Rs")
    sd_via_R = load(sfile,"sd_via_R")
    sd_ηs_R = load(sfile,"sd_ηs_R")
    sd_ωs_R = load(sfile,"sd_ωs_R")
    sd_kc_R = load(sfile,"sd_kc_R")
    sd_KS_R = load(sfile,"sd_KS_R")
    sd_kr_R = load(sfile,"sd_kr_R")
    sd_ηs = load(sfile,"sd_ηs")
    sd_via_η = load(sfile,"sd_via_η")
    sd_ωs = load(sfile,"sd_ωs")
    sd_via_ω = load(sfile,"sd_via_ω")
    sd_fr_ΔG = load(sfile,"sd_fr_ΔG")
    sd_via_a = load(sfile,"sd_via_a")
    sd_via_ϕR = load(sfile,"sd_via_ϕR")
    sd_kcs = load(sfile,"sd_kcs")
    sd_KSs = load(sfile,"sd_KSs")
    sd_krs = load(sfile,"sd_krs")
    sd_av_steps = load(sfile,"sd_av_steps")
    # Preallocate standard errors
    se_Rs = zeros(size(sd_Rs))
    se_via_R = zeros(size(sd_via_R))
    se_ηs_R = zeros(size(sd_ηs_R))
    se_ωs_R = zeros(size(sd_ωs_R))
    se_kc_R = zeros(size(sd_kc_R))
    se_KS_R = zeros(size(sd_KS_R))
    se_kr_R = zeros(size(sd_kr_R))
    # Calculate standard errors from this
    se_svt = sd_svt./sqrt.(no_sims)
    se_tsvt = sd_tsvt./sqrt.(no_sims)
    se_sbs = sd_sbs./sqrt.(no_sims)
    se_ηs = sd_ηs./sqrt.(no_sims)
    se_ωs = sd_ωs./sqrt.(no_sims)
    # Calculation (slightly) different in the viable case
    se_via_η = sd_via_η./sqrt.(no_via)
    se_via_ω = sd_via_ω./sqrt.(no_via)
    se_fr_ΔG = sd_fr_ΔG./sqrt.(no_via)
    se_via_a = sd_via_a./sqrt.(no_via)
    se_via_ϕR = sd_via_ϕR./sqrt.(no_via)
    se_kcs = sd_kcs./sqrt.(no_via)
    se_KSs = sd_KSs./sqrt.(no_via)
    se_krs = sd_krs./sqrt.(no_via)
    se_av_steps = sd_av_steps./sqrt.(no_via)
    for i = 1:size(sd_Rs,1)
        se_Rs[i,:] = sd_Rs[i,:]./sqrt.(no_sims)
        se_via_R[i,:] = sd_via_R[i,:]./sqrt.(no_via)
        se_ηs_R[i,:] = sd_ηs_R[i,:]./sqrt.(no_rs[i,:])
        se_ωs_R[i,:] = sd_ωs_R[i,:]./sqrt.(no_rs[i,:])
        se_kc_R[i,:] = sd_kc_R[i,:]./sqrt.(no_rs[i,:])
        se_KS_R[i,:] = sd_KS_R[i,:]./sqrt.(no_rs[i,:])
        se_kr_R[i,:] = sd_kr_R[i,:]./sqrt.(no_rs[i,:])
    end
    # Check if directory exists and if not make it
    if ~isdir("Output/Plotsd=$(d)u=$(μrange)")
        mkdir("Output/Plotsd=$(d)u=$(μrange)")
    end
    if sim_type == 1
        tl = "high free-energy low loss"
    elseif sim_type == 2
        tl = "low free-energy low loss"
    elseif sim_type == 3
        tl = "high free-energy high loss"
    elseif sim_type == 4
        tl = "low free-energy high loss"
    end
    # Setup plotting
    pyplot(dpi=200,title=tl)
    plot(xlabel="Time (s)",ylabel="Number of surviving strains",xlim=(-Inf,5e7))
    plot!(times,mn_svt,ribbon=se_svt,label="")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvSurvTime.png")
    plot(xlabel="Time (s)",ylabel="Number of viable strains",xlim=(-Inf,5e7))
    plot!(times,mn_tsvt,ribbon=se_tsvt,label="")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvViaTime.png")
    plot(xlabel="Time (s)",ylabel="Number of diversified substrates",xlim=(-Inf,5e7))
    plot!(times,mn_sbs,ribbon=se_sbs,label="")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvSubTime.png")
    plot(xlabel="Time (s)",ylabel="Number of strains",xlim=(-Inf,5e7))
    plot!(times,mn_Rs[1,:],ribbon=se_Rs[1,:],label="R=1")
    plot!(times,mn_Rs[3,:],ribbon=se_Rs[3,:],label="R=3")
    plot!(times,mn_Rs[5,:],ribbon=se_Rs[5,:],label="R=5")
    plot!(times,mn_Rs[7,:],ribbon=se_Rs[7,:],label="R=7")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvReacsTime.png")
    plot(xlabel="Time (s)",ylabel="Number of strains",xlim=(-Inf,5e7))
    plot!(times,mn_via_R[1,:],ribbon=se_via_R[1,:],label="R=1")
    plot!(times,mn_via_R[3,:],ribbon=se_via_R[3,:],label="R=3")
    plot!(times,mn_via_R[5,:],ribbon=se_via_R[5,:],label="R=5")
    plot!(times,mn_via_R[7,:],ribbon=se_via_R[7,:],label="R=7")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvViaReacsTime.png")
    plot(xlabel="Time (s)",ylabel="Average eta value",xlim=(-Inf,5e7))
    plot!(times,mn_ηs_R[1,:],ribbon=se_ηs_R[1,:],label="R=1")
    plot!(times,mn_ηs_R[3,:],ribbon=se_ηs_R[3,:],label="R=3")
    plot!(times,mn_ηs_R[5,:],ribbon=se_ηs_R[5,:],label="R=5")
    plot!(times,mn_ηs_R[7,:],ribbon=se_ηs_R[7,:],label="R=7")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvEtaperReacTime.png")
    plot(xlabel="Time (s)",ylabel="Average omega value",xlim=(-Inf,5e7))
    plot!(times,mn_ωs_R[1,:],ribbon=se_ωs_R[1,:],label="R=1")
    plot!(times,mn_ωs_R[3,:],ribbon=se_ωs_R[3,:],label="R=3")
    plot!(times,mn_ωs_R[5,:],ribbon=se_ωs_R[5,:],label="R=5")
    plot!(times,mn_ωs_R[7,:],ribbon=se_ωs_R[7,:],label="R=7")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvOmegaperReacTime.png")
    plot(xlabel="Time (s)",ylabel="Average forward rate constant",xlim=(-Inf,5e7))
    plot!(times,mn_kc_R[1,:],ribbon=se_kc_R[1,:],label="R=1")
    plot!(times,mn_kc_R[3,:],ribbon=se_kc_R[3,:],label="R=3")
    plot!(times,mn_kc_R[5,:],ribbon=se_kc_R[5,:],label="R=5")
    plot!(times,mn_kc_R[7,:],ribbon=se_kc_R[7,:],label="R=7")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvkcperReacTime.png")
    plot(xlabel="Time (s)",ylabel="Average half saturation constant",xlim=(-Inf,5e7))
    plot!(times,mn_KS_R[1,:],ribbon=se_KS_R[1,:],label="R=1")
    plot!(times,mn_KS_R[3,:],ribbon=se_KS_R[3,:],label="R=3")
    plot!(times,mn_KS_R[5,:],ribbon=se_KS_R[5,:],label="R=5")
    plot!(times,mn_KS_R[7,:],ribbon=se_KS_R[7,:],label="R=7")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvKSperReacTime.png")
    plot(xlabel="Time (s)",ylabel="Average reversibility factor",xlim=(-Inf,5e7))
    plot!(times,mn_kr_R[1,:],ribbon=se_kr_R[1,:],label="R=1")
    plot!(times,mn_kr_R[3,:],ribbon=se_kr_R[3,:],label="R=3")
    plot!(times,mn_kr_R[5,:],ribbon=se_kr_R[5,:],label="R=5")
    plot!(times,mn_kr_R[7,:],ribbon=se_kr_R[7,:],label="R=7")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvkrperReacTime.png")
    plot(xlabel="Time (s)",ylabel="Average eta value",xlim=(-Inf,5e7))
    plot!(times,mn_ηs,ribbon=se_ηs,label="")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvEtaTime.png")
    plot(xlabel="Time (s)",ylabel="Average eta value",xlim=(-Inf,5e7))
    plot!(times,mn_via_η,ribbon=se_via_η,label="")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvViaEtaTime.png")
    plot(xlabel="Time (s)",ylabel="Average omega value",xlim=(-Inf,5e7))
    plot!(times,mn_ωs,ribbon=se_ωs,label="")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvOmegaTime.png")
    plot(xlabel="Time (s)",ylabel="Average omega value",xlim=(-Inf,5e7))
    plot!(times,mn_via_ω,ribbon=se_via_ω,label="")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvViaOmegaTime.png")
    plot(xlabel="Time (s)",ylabel="Fraction transduced",xlim=(-Inf,5e7))
    plot!(times,mn_fr_ΔG,ribbon=se_fr_ΔG,label="")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvFracTransTime.png")
    plot(xlabel="Time (s)",ylabel="Average forward rate",xlim=(-Inf,5e7))
    plot!(times,mn_kcs,ribbon=se_kcs,label="")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvkcTime.png")
    plot(xlabel="Time (s)",ylabel="Average saturation constant",xlim=(-Inf,5e7))
    plot!(times,mn_KSs,ribbon=se_KSs,label="")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvKSTime.png")
    plot(xlabel="Time (s)",ylabel="Average reversibility factor",xlim=(-Inf,5e7))
    plot!(times,mn_krs,ribbon=se_krs,label="")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvkrTime.png")
    plot(xlabel="Time (s)",ylabel="Average steps in metabolite hierachy",xlim=(-Inf,5e7))
    plot!(times,mn_av_steps,ribbon=se_av_steps,label="")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvStepsTime.png")
    plot(xlabel="Time (s)",ylabel="Average ATP concentration",xlim=(-Inf,5e7))
    plot!(times,mn_via_a,ribbon=se_via_a,label="")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvATPTime.png")
    plot(xlabel="Time (s)",ylabel="Average Ribosome fraction",xlim=(-Inf,5e7))
    plot!(times,mn_via_ϕR,ribbon=se_via_ϕR,label="")
    savefig("Output/Plotsd=$(d)u=$(μrange)/AvFracTime.png")
    return(nothing)
end

# Function to load in and plot the average interaction strengths
function plot_av_ints()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 1
        error("insufficent inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    ims = 0
    # Check that all arguments can be converted to integers
    try
        ims = parse(Int64,ARGS[1])
    catch e
        error("need to provide an integer")
    end
    println("Compiled")
    # Load in hardcoded simulation parameters
    Np, Nt, M, d = sim_paras()
    # Find file name to load in
    sfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/IntStats$(ims)Ims.jld"
    # Check it actually exists
    if ~isfile(sfile)
        error("missing stats file for $(ims) immigrations simulations")
    end
    # Now load out the times, and number of trajectories
    times = load(sfile,"times")
    no_sims = load(sfile,"no_sims")
    no_via = load(sfile,"no_via")
    # Load in averages
    mn_svt = load(sfile,"mn_svt")
    mn_tsvt = load(sfile,"mn_tsvt")
    mn_no_comp = load(sfile,"mn_no_comp")
    mn_no_facl = load(sfile,"mn_no_facl")
    mn_no_selff = load(sfile,"mn_no_selff")
    mn_no_selfc = load(sfile,"mn_no_selfc")
    mn_via_no_comp = load(sfile,"mn_via_no_comp")
    mn_via_no_facl = load(sfile,"mn_via_no_facl")
    mn_via_no_selff = load(sfile,"mn_via_no_selff")
    mn_via_no_selfc = load(sfile,"mn_via_no_selfc")
    mn_st_comp = load(sfile,"mn_st_comp")
    mn_st_facl = load(sfile,"mn_st_facl")
    mn_st_selfc = load(sfile,"mn_st_selfc")
    mn_st_selff = load(sfile,"mn_st_selff")
    # Load in standard deviations
    sd_svt = load(sfile,"sd_svt")
    sd_tsvt = load(sfile,"sd_tsvt")
    sd_no_comp = load(sfile,"sd_no_comp")
    sd_no_facl = load(sfile,"sd_no_facl")
    sd_no_selff = load(sfile,"sd_no_selff")
    sd_no_selfc = load(sfile,"sd_no_selfc")
    sd_via_no_comp = load(sfile,"sd_via_no_comp")
    sd_via_no_facl = load(sfile,"sd_via_no_facl")
    sd_via_no_selff = load(sfile,"sd_via_no_selff")
    sd_via_no_selfc = load(sfile,"sd_via_no_selfc")
    sd_st_comp = load(sfile,"sd_st_comp")
    sd_st_facl = load(sfile,"sd_st_facl")
    sd_st_selfc = load(sfile,"sd_st_selfc")
    sd_st_selff = load(sfile,"sd_st_selff")
    # Calculate standard errors from this
    se_svt = sd_svt./sqrt.(no_sims)
    se_tsvt = sd_tsvt./sqrt.(no_sims)
    se_no_comp = sd_no_comp./sqrt.(no_sims)
    se_no_facl = sd_no_facl./sqrt.(no_sims)
    se_no_selff = sd_no_selff./sqrt.(no_sims)
    se_no_selfc = sd_no_selfc./sqrt.(no_sims)
    se_st_comp = sd_st_comp./sqrt.(no_sims)
    se_st_facl = sd_st_facl./sqrt.(no_sims)
    se_st_selff = sd_st_selff./sqrt.(no_sims)
    se_st_selfc = sd_st_selfc./sqrt.(no_sims)
    # Calculation (slightly) different in the viable case
    se_via_no_comp = sd_via_no_comp./sqrt.(no_via)
    se_via_no_facl = sd_via_no_facl./sqrt.(no_via)
    se_via_no_selff = sd_via_no_selff./sqrt.(no_via)
    se_via_no_selfc = sd_via_no_selfc./sqrt.(no_via)
    # Setup plotting
    # ALL QUITE MUCKY HERE AT THE MOMENT, CLEAN THESE UP BEFORE I PRESENT THEM
    pyplot(dpi=200)
    plot(xlabel="Time (s)",ylabel="Number of competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_no_comp./(mn_no_selfc.^2),label="")
    savefig("Output/AvCompTime.png")
    plot(xlabel="Time (s)",ylabel="Number of faciliation interactions",xlim=(-Inf,5e7))
    plot!(times,mn_no_facl./(mn_no_selfc.^2),label="")
    savefig("Output/AvFaclTime.png")
    plot(xlabel="Time (s)",ylabel="Number of (self) competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_no_selfc./mn_no_selfc,label="")
    savefig("Output/AvSelfCompTime.png")
    plot(xlabel="Time (s)",ylabel="Number of (self) faciliation interactions",xlim=(-Inf,5e7))
    plot!(times,mn_no_selff./mn_no_selfc,label="")
    savefig("Output/AvSelfFaclTime.png")
    plot(xlabel="Time (s)",ylabel="Number of viable competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_via_no_comp./(mn_via_no_selfc.^2),label="")
    savefig("Output/AvViaCompTime.png")
    plot(xlabel="Time (s)",ylabel="Number of viable faciliation interactions",xlim=(-Inf,5e7))
    plot!(times,mn_via_no_facl./(mn_via_no_selfc.^2),label="")
    savefig("Output/AvViaFaclTime.png")
    plot(xlabel="Time (s)",ylabel="Number of viable (self) competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_via_no_selfc./mn_via_no_selfc,label="")
    savefig("Output/AvViaSelfCompTime.png")
    plot(xlabel="Time (s)",ylabel="Number of viable (self) competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_via_no_selff./mn_via_no_selfc,label="")
    savefig("Output/AvViaSelfFaclTime.png")
    plot(xlabel="Time (s)",ylabel="Strength of competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_st_comp,ribbon=se_st_comp,label="")
    savefig("Output/AvStCompTime.png")
    plot(xlabel="Time (s)",ylabel="Strength of faciliation interactions",xlim=(-Inf,5e7))
    plot!(times,mn_st_facl,ribbon=se_st_facl,label="")
    savefig("Output/AvStFaclTime.png")
    plot(xlabel="Time (s)",ylabel="Strength of (self) competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_st_selfc,ribbon=se_st_selfc,label="")
    savefig("Output/AvStSelfCompTime.png")
    plot(xlabel="Time (s)",ylabel="Strength of (self) competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_st_selff,ribbon=se_st_selff,label="")
    savefig("Output/AvStSelfFaclTime.png")
    plot(xlabel="Time (s)",ylabel="Mean strength of competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_st_comp./mn_no_comp,label="")
    savefig("Output/MnStCompTime.png")
    plot(xlabel="Time (s)",ylabel="Mean strength of faciliation interactions",xlim=(-Inf,5e7))
    plot!(times,mn_st_facl./mn_no_facl,label="")
    savefig("Output/MnStFaclTime.png")
    plot(xlabel="Time (s)",ylabel="Mean strength of (self) competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_st_selfc./mn_no_selfc,label="")
    savefig("Output/MnStSelfCompTime.png")
    plot(xlabel="Time (s)",ylabel="Mean strength of (self) competition interactions",xlim=(-Inf,5e7))
    plot!(times,mn_st_selff./mn_no_selff,label="")
    savefig("Output/MnStSelfFaclTime.png")
    return(nothing)
end

# function to plot the averages for just one run
function plot_run_averages()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 2
        error("insufficent inputs provided")
    end
    # Preallocate the variables I want to extract from the input
    ims = 0
    rN = 0
    # Check that all arguments can be converted to integers
    try
        ims = parse(Int64,ARGS[1])
        rN = parse(Int64,ARGS[2])
    catch e
            error("need to provide two integers")
    end
    println("Compiled")
    # Load in hardcoded simulation parameters
    Np, Nt, M, d = sim_paras()
    # Read in appropriate files
    pfile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/Paras$(ims)Ims.jld"
    if ~isfile(pfile)
        error("$(ims) immigrations run $(rN) is missing a parameter file")
    end
    ofile = "Output/$(Np)Pools$(M)Metabolites$(Nt)Speciesd=$(d)/AvRun$(rN)Data$(ims)Ims.jld"
    if ~isfile(ofile)
        error("$(ims) immigrations run $(rN) is missing an output file")
    end
    # Read in relevant data
    ps = load(pfile,"ps")
    T = load(ofile,"T")
    svt = load(ofile,"svt")
    tsvt = load(ofile,"tsvt")
    ηs = load(ofile,"ηs")
    # Setup plotting
    pyplot(dpi=200)
    # Plot this data
    plot(T,svt,label="",xlabel="Time (s)",ylabel="Number of survivors")
    savefig("Output/SvTime.png")
    plot(T,ηs,label="",xlabel="Time (s)",ylabel="eta")
    savefig("Output/EtaTime.png")
    plot(T,tsvt,label="",xlabel="Time (s)",ylabel="Number of viable strains")
    savefig("Output/ViaTime.png")
    return(nothing)
end

# @time plot_run_averages()
# @time plot_aves()
@time plot_traj()
# @time plot_av_ints()
