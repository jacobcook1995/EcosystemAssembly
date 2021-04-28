# Script to construct figure 5
using Assembly
using Plots
using JLD
using LaTeXStrings
using StatsBase
using Plots.PlotMeasures
using LsqFit
import PyPlot

# Make figure plots for the interactions
function SI_ints(Rl::Int64,Ru::Int64,syns::Array{Bool,1},rps::Int64,Ni::Int64,en::String)
    println("Compiled!")
    # Preallocate memory to store number of interactions
    ins1 = zeros(rps,length(syns))
    ins2 = zeros(rps,length(syns))
    ins3 = zeros(rps,length(syns))
    ins4 = zeros(rps,length(syns))
    # Preallocate memory to store mean interaction strengths
    mean1 = fill(NaN,rps,length(syns))
    mean2 = fill(NaN,rps,length(syns))
    mean3 = fill(NaN,rps,length(syns))
    mean4 = fill(NaN,rps,length(syns))
    # Preallocate vectors to store all interaction strengths
    sts1 = Array{Array{Float64,1},1}(undef,2)
    sts2 = Array{Array{Float64,1},1}(undef,2)
    sts3 = Array{Array{Float64,1},1}(undef,2)
    sts4 = Array{Array{Float64,1},1}(undef,2)
    # Make array of bools to check if data has already been saved
    fll = fill(false,4,2)
    # Loop over parameter sets
    for i = 1:rps
        for j = 1:length(syns)
            # Read in relevant files
            ifile = "Data/$(Rl)-$(Ru)$(syns[j])$(Ni)$(en)/IntsReacs$(Rl)-$(Ru)Syn$(syns[j])Run$(i)Ns$(Ni).jld"
            if ~isfile(ifile)
                error("run $(i) is missing an interaction file")
            end
            # Just need to load ints and their strengths
            ints = load(ifile,"ints")
            in_str = load(ifile,"in_str")
            # Store number of each interaction type
            ins1[i,j] = count(x->x==1,ints)
            ins2[i,j] = count(x->x==2,ints)
            ins3[i,j] = count(x->x==3,ints)
            ins4[i,j] = count(x->x==4,ints)
            # Find corresponding positions
            pos1 = findall(x->x==1,ints)
            pos2 = findall(x->x==2,ints)
            pos3 = findall(x->x==3,ints)
            pos4 = findall(x->x==4,ints)
            # Use to find vector of interaction strengths
            str1 = in_str[pos1]
            str2 = in_str[pos2]
            str3 = in_str[pos3]
            str4 = in_str[pos4]
            # Check if there are any interactions of that type
            if length(str1) >= 1
                # Find mean strength of each interaction
                mean1[i,j] = mean(str1)
                # If vector doesn't already exist create it
                if fll[1,j] == false
                    sts1[j] = str1
                    # Update counter
                    fll[1,j] = true
                else
                    sts1[j] = cat(sts1[j],str1,dims=1)
                end
            end
            # Do for each interaction type
            if length(str2) >= 1
                mean2[i,j] = mean(str2)
                # If vector doesn't already exist create it
                if fll[2,j] == false
                    sts2[j] = str2
                    # Update counter
                    fll[2,j] = true
                else
                    sts2[j] = cat(sts2[j],str2,dims=1)
                end
            end
            if length(str3) >= 1
                mean3[i,j] = mean(str3)
                # If vector doesn't already exist create it
                if fll[3,j] == false
                    sts3[j] = str3
                    # Update counter
                    fll[3,j] = true
                else
                    sts3[j] = cat(sts3[j],str3,dims=1)
                end
            end
            if length(str4) >= 1
                mean4[i,j] = mean(str4)
                # If vector doesn't already exist create it
                if fll[4,j] == false
                    sts4[j] = str4
                    # Update counter
                    fll[4,j] = true
                else
                    sts4[j] = cat(sts4[j],str4,dims=1)
                end
            end
        end
    end
    # Find percentage of each interactions type
    cvsf = 100.0*(ins1)./(ins1.+ins2)
    pvsn = 100.0*(ins2.+ins3)./(ins1.+ins2.+ins3.+ins4)
    tvnt = 100.0*(ins3.+ins4)./(ins1.+ins2.+ins3.+ins4)
    # Setup plotting
    pyplot()
    theme(:wong2,dpi=300)
    # Preallocate array to store plots
    p = Array{Plots.Plot,2}(undef,2,length(syns))
    # Loop over syntrophy conditions
    for j = 1:length(syns)
        # Make plot title
        tl = ""
        if syns[j] == true
            tl = "Reversible kinetics"
        else
            tl = "Michaelis–Menten kinetics"
        end
        # Make range of ticks to label
        rgn = collect(-15:3:0)
        ergn = fill("",length(rgn))
        # Find corresponding exponentials
        for i = 1:length(rgn)
            ergn[i] = L"10^{%$(rgn[i])}"
        end
        # Make bins for proportions
        sbins = range(-15.0,stop=0.0,length=52)
        # Now plot all strengths
        p[1,j] = plot(title=tl,ylabel="Number of interactions",legend=:topleft)
        plot!(p[1,j],xticks=(rgn,ergn),legendfontsize=10,tickfontsize=10,guidefontsize=12)
        # Add xlabel for reversible case
        if syns[j] == true
            plot!(p[1,j],legend=false,xlabel="Interaction strength")
        end
        histogram!(p[1,j],log10.(sts1[j]),fillalpha=0.75,label="Competiton",bins=sbins)
        histogram!(p[1,j],log10.(sts2[j]),fillalpha=0.75,label="Facilitation",bins=sbins)
        histogram!(p[1,j],log10.(sts4[j]),fillalpha=0.75,label="Pollution",bins=sbins)
        histogram!(p[1,j],log10.(sts3[j]),fillalpha=0.75,label="Syntrophy",bins=sbins)
        # Choose which letter to annotate
        if syns[j] == false
            # Add annotation
            px, py = annpos([-15.0,0.0],[0.0,4300.0],0.15,0.05)
            annotate!(p[1,j],px,py,text("A",17,:black))
        else
            # Add annotation
            px, py = annpos([-15.0,0.0],[0.0,4100.0],0.15,0.05)
            annotate!(p[1,j],px,py,text("C",17,:black))
        end
        # Fit pollution histogram
        hp = fit(Histogram,log10.(sts4[j]),sbins,closed=:right)
        # Find maximum height of this pollution histogram
        hmaxp = maximum(hp.weights)
        # Find index of this maximum
        mind = findfirst(x->x==hmaxp,hp.weights)
        # Find first index after this that's less than half the value
        dind = findfirst(x->x<=0.5*hmaxp,hp.weights[(mind+1):end])
        # Extract bins from histogram
        rn = collect(sbins)
        # Next need to find x posistion of this half maximum point
        xpos = (rn[mind+dind]+rn[mind+dind-1])/2
        # Put in a vertical line here
        vline!(p[1,j],[xpos],color=:red,label="")
        savefig(p[1,j],"Output/SI/AllIntStrength$(Rl)-$(Ru)$(syns[j])$(Ni)$(en).png")
        # Move onto plotting the simplex
        p[2,j] = plot(grid=false,showaxis=false,xlim=(-0.325,1.325),ylim=(-0.1,1.0))
        plot!(p[2,j],legendfontsize=10,tickfontsize=10,guidefontsize=12)
        # Plot the triangle over this
        plot!(p[2,j],[0.0;0.5],[0.0;0.8660],color=:black,label=false)
        plot!(p[2,j],[0.0;1.0],[0.0;0.0],color=:black,label=false)
        plot!(p[2,j],[1.0;0.5],[0.0;0.8660],color=:black,label=false)
        # Bottom left competition
        annotate!(p[2,j],0.0,-0.04,text("Competition",12,:black))
        # Bottom right thermodynamic
        annotate!(p[2,j],1.0,-0.04,text("Facilitation",12,:black))
        # Top facilitation
        annotate!(p[2,j],0.5,0.9,text("Thermodynamic",12,:black))
        # Add xlabel for reversible case
        if syns[j] == true
            annotate!(p[2,j],0.5,-0.21,text("Proportion of interactions",12,:black))
        end
        # Group by interaction type
        a = ins1[:,j] # Competiton
        b = ins2[:,j] # Facilitation
        c = ins3[:,j] .+ ins4[:,j] # Thermodynamic
        # Find x and y coordinates for each point
        x = 0.5*(2 .*b.+c)./(a.+b.+c)
        y = (sqrt(3)/2)*(c)./(a.+b.+c)
        scatter!(p[2,j],x,y,label="")
        # Choose which letter to annotate
        if syns[j] == false
            # Add annotation
            annotate!(p[2,j],-0.15,1.0,text("B",17,:black))
        else
            # Add annotation
            annotate!(p[2,j],-0.15,1.0,text("D",17,:black))
        end
        # Then save the figure
        savefig(p[2,j],"Output/SI/Simplex$(Rl)-$(Ru)$(syns[j])$(Ni)$(en).png")
    end
    # Combine all graphs and save
    pt = plot(p[1,2],p[2,2],p[1,1],p[2,1],layout=4,size=(1200,800),margin=5.0mm)
    # Save as high energy supply case
    savefig(pt,"Output/SI/HighInts.png")
    return(nothing)
end

# function to plot strain diversity against substrate diversity
function SvvsDv(Nr::Int64,Ni::Int64,Rls::Array{Int64,1},Rus::Array{Int64,1},syns::Array{Bool,1},
                ens::Array{String,1})
    # Find number of parameter sets
    Ns = length(Rls)
    # Container to store number of survivors
    svs = zeros(Int64,Ns,Nr)
    # Container to store metabolite diversity
    mbs = zeros(Int64,Ns,Nr)
    # Loop over parameter sets
    for i = 1:length(Rls)
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
            for k = 2:ps.M
                if inf_out[ps.N+k] > 0.0
                    # Increment counter
                    cm += 1
                end
            end
            # Save results to vector
            mbs[i,j] = cm
        end
    end
    # Setup plotting
    pyplot()
    theme(:wong2,dpi=200)
    wongc = wong2_palette()
    # Make scatter plot of substrate diversification against strain diversity
    plot(title="Survivors per substrate",ylabel="Number of survivors",xlabel="Number of substrates")
    for i = 1:Ns
        scatter!(mbs[i,:],svs[i,:],color=wongc[1],label="")
    end
    # Plot maximum line
    plot!(1:24,1:24,color=:red,label="")
    savefig("Output/SI/SubstrateDiversity.png")
    return(nothing)
end

# Function to make growth laws plot
function growth_laws()
    println("Successfully compiled.")
    # Simple test data set
    ai = 5.0 # initial energy level
    Ni = 100.0 # initial population
    Si = 1.0
    Pi = 0.0
    ϕi = 0.1
    # Choose simulation time
    Tmax = 5000000.0
    # Set this as a middling value of ΔG
    ΔG = -3e6 # Relatively small Gibbs free energy change
    # Max Elongation rate also taken from Bremer (1996), convert from minutes to seconds
    γm = 1260.0/60.0 # Change this one
    # Make vector of γm values
    γs = [γm/100,γm/50,γm/20,γm/10,γm/7.5,γm/5,γm/4,γm/3,γm/2,γm/1.5,γm]
    # Make vector to store final growth rates and fractions
    λ1 = zeros(length(γs))
    ϕ1 = zeros(length(γs))
    a1 = zeros(length(γs))
    # Setup plotting options
    pyplot(dpi=200)
    pR = L"\phi_R"
    p1 = plot(xaxis="Growth rate, λ",yaxis="Ribosome fraction, $(pR)")
    for i = 1:length(γs)
        ps = initialise_prot_gl(γs[i],ΔG)
        C, T = prot_simulate_fix(ps,Tmax,ai,Ni,Si,Pi,ϕi)
        # Now calculate growth rates and proteome fractions
        λa = zeros(length(T))
        ϕR = zeros(length(T))
        for j = 1:length(T)
            ϕR[j] = ϕ_R(C[j,2],ps)
            λa[j] = λs(C[j,2],ϕR[j],ps)
        end
        # Save final λ and ϕR values
        λ1[i] = λa[end]
        ϕ1[i] = ϕR[end]
        a1[i] = C[end,2]
    end
    # Now make set of ΔG values
    ΔGs = [ΔG,ΔG/1.5,ΔG/2,ΔG/3,ΔG/4,ΔG/5,ΔG/7.5,ΔG/10,ΔG/20,ΔG/50,ΔG/100]
    # Make vector to store final growth rates and fractions
    λ2 = zeros(length(ΔGs))
    ϕ2 = zeros(length(ΔGs))
    a2 = zeros(length(ΔGs))
    for i = 1:length(ΔGs)
        ps = initialise_prot_gl(γm,ΔGs[i])
        C, T = prot_simulate_fix(ps,Tmax,ai,Ni,Si,Pi,ϕi)
        # Now calculate growth rates and proteome fractions
        λa = zeros(length(T))
        ϕR = zeros(length(T))
        for j = 1:length(T)
            ϕR[j] = ϕ_R(C[j,2],ps)
            λa[j] = λs(C[j,2],ϕR[j],ps)
        end
        # Save final λ and ϕR values
        λ2[i] = λa[end]
        ϕ2[i] = ϕR[end]
        a2[i] = C[end,2]
    end
    # Now want to do a least squares fit for both sets of data
    @. model(x, p) = p[1]*x + p[2]
    p0 = [0.5,0.5]
    fit1 = curve_fit(model,λ1[4:end],ϕ1[4:end],p0)
    pr1 = coef(fit1)
    fit2 = curve_fit(model,λ2[1:end-3],ϕ2[1:end-3],p0)
    pr2 = coef(fit2)
    # plot both lines on the graph
    λ1s = [0.0;λ1]
    λ2s = [0.0;λ2]
    plot!(p1,λ1s,pr1[1]*λ1s.+pr1[2],label="",color=:blue)
    plot!(p1,λ2s,pr2[1]*λ2s.+pr2[2],label="",color=:red)
    # Add arrows indicating direction of change
    l = 1e-4 # Way too large
    quiver!(p1,[λ1[end-2]],[ϕ1[end-2]+0.03],quiver=([-l],[-pr1[1]*l]),color=:blue)
    quiver!(p1,[λ2[6]],[ϕ2[6]-0.03],quiver=([l],[pr2[1]*l]),color=:red)
    # Position is where the annotation centres are
    pos1x = λ1[end-2] - l/2
    pos1y = ϕ1[end-2] + 0.05 - pr1[1]*l/2
    pos2x = λ2[6] + l/2
    pos2y = ϕ2[6] - 0.05 + pr2[1]*l/2
    # Calculate rotations in degrees
    r1 = -12.5 # Needs to be -ve
    r2 = 20 # Needs to be +ve
    # Then add the annotations
    annotate!(p1,pos1x,pos1y,text("Translational inhibition",8,color=:blue,rotation=r1))
    annotate!(p1,pos2x,pos2y,text("Nutrient quality",8,color=:red,rotation=r2))
    # Plot final values
    scatter!(p1,λ1,ϕ1,label="",color=:lightblue,markersize=5)
    scatter!(p1,λ2,ϕ2,label="",color=:orange,markersize=5)
    # Finally save the graph
    savefig(p1,"Output/SI/GrowthLaws.png")
    return(nothing)
end

# Run both high interactions, survivors vs diversity, and growth law plots
@time SI_ints(1,5,[true,false],250,250,"h")
@time SvvsDv(250,250,[1,1,1,1],[5,5,5,5],[false,true,false,true],["l","l","h","h"])
@time growth_laws()
