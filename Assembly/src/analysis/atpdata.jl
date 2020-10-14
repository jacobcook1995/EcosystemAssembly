# Script to read in organise and plot the ATP data
using Assembly
using CSV
using DataFrames
using Plots
using Statistics
using StatsBase
using DifferentialEquations
import PyPlot

# Writing a function to perform a hampel filter on a data set
function hampel_filter(x::Array{Float64,1},k::Int64)
    # Preallocate vector of points to keep
    kp = fill(true,length(x))
    # First and last point are preserved
    for i = 2:(length(x)-1)
        # Find start point of the sample
        if i <= k
            st = 1
        else
            st = i - k
        end
        # Find end point of the sample
        if i + k >= length(x)
            ed = length(x)
        else
            ed = i + k
        end
        # Now make sample
        sx = x[st:ed]
        # Find median of the sample
        mx = median(sx)
        # Find median absolute deviation of the sample about the median
        mdx = mad(x;center=mx,normalize=false)
        # Convert into standard deviation
        sdx = 1.4826*mdx
        # Check if point lies within 3 standard deviations of sample median
        if x[i] > mx + 3*sdx || x[i] < mx - 3*sdx
            kp[i] = false
        end
    end
    return(kp)
end

# function to evaluate energy concentration with time
function dadt!(dx::Array{Float64,1},x::Array{Float64,1},m::MicrobeP,S::Float64,P::Float64,t::Float64)
    # Need to find E
    E = Eα(x[2],m,1)
    # Assume far from equilbrium
    θs = 0.0
    # Find rate of energy aquisition
    J = m.η[1]*qs(m,S,P,E,θs)
    # Find elongation rate
    γ = γs(x[1],m)
    # Calculate growth rate
    λ = λs(x[1],x[2],m)
    # Calculate total change in energy concentration
    dx[1] = J - (m.ρ*m.MC+x[1])*λ
    # Calculate optimal ribosome fraction
    ϕR_opt = ϕ_R(x[1],m)
    # Find time delay
    τ = m.fd/λ
    # Also update ribosome fraction
    dx[2] = (ϕR_opt - x[2])/τ
    return(dx)
end

function eval_a(m::MicrobeP,S::Float64,P::Float64,ϕR0::Float64,Tmax::Float64,tsave::Array{Float64,1})
    # substitute constants into expression
    dfdt!(dx,x,m,t) = dadt!(dx,x,m,S,P,t)
    # Find time span for this step
    tspan = (0,Tmax)
    x0 = [10.0;ϕR0]
    # Then setup and solve the problem
    prob = ODEProblem(dfdt!,x0,tspan,m)
    # Still generates problems, not sure if I have to change a solver option or what
    sol = DifferentialEquations.solve(prob,saveat=tsave)
    # The biggest issue is how to find this at particular time points
    return(sol',sol.t)
end

# function to find chi_squared difference between data and simulation
function chi_square(m::MicrobeP,S::Float64,P::Float64,ϕR0::Float64,Tmax::Float64,tsave::Array{Float64,1},
                    indm::Int64,admA::Array{Float64,1},adsdA::Array{Float64})
    # Run function to evaluate solution
    C, T = eval_a(m,S,P,ϕR0,Tmax,tsave)
    # Calculate means squared error here now
    χ2 = 0
    # Set range and overwrite if necessary
    rang = 1:length(T)
    if indm == 0
        rang = 2:length(T)
    end
    for i = rang
        χ2 += ((C[i,1]-admA[indm+(i-1)])^2)/((adsdA[indm+(i-1)])^2)
    end
    return(χ2)
end

# function to use steppest descent to find optimal paramters to fit the emperical curves
function prot_fit(times::Array{Float64,1},pt::Float64,mA::Array{Float64,1},sdA::Array{Float64,1})
    # Change units to molecules per cell
    admA = mA*6.02e23/1e9
    adsdA = sdA*6.02e23/1e9
    # Choose initial ϕR value (not treating as a parameter)
    ϕR0 = 0.05
    # Choose concentration of glucose in blood as substrate concentration
    S = 5.5e-3
    # Product fixed at low value as not interested in inhibition here
    P = S/100.0
    # Use parameters from literature (justifuied in make_full.jl)
    MC = 10^8
    γm = 1260.0/60.0
    n = zeros(Int64,3)
    n[1] = 7459
    n[2:3] .= 300
    ρ = 29.0
    Pb = 0.7
    ϕH = 0.45
    fd = log(100)/log(2)
    # death rate isn't actually going to be used
    d = 6.0e-5
    # Only considering one reaction for simplicity
    R = 1
    Reacs = [1]
    kcs = [10.0]
    KSs = [(1/4)*5.5e-3]
    krs = [10.0]
    ϕP = [1.0]
    # No explict reactions so η is a parameter with less consequence
    η = [1.0]
    # Calculate time span in seconds
    Tmax = (maximum(times)-pt)*60.0
    # Find index of value corresponding to peak
    _, indm = findmax(mA)
    # Find value one step before the peak
    indm -= 1
    # Select time points to save at
    tsave = (times .- pt)*60.0
    # As an initial test I'm going to try these parameter choices
    Kγ = 5e8
    KΩ = 1e9
    m = make_MicrobeP(MC,γm,ρ,Kγ,Pb,d,ϕH,KΩ,fd,R,Reacs,η,kcs,KSs,krs,n,ϕP)
    # Use function to find χ2
    χ2 = chi_square(m,S,P,ϕR0,Tmax,tsave,indm,admA,adsdA)
    # Set up while loop
    stab = false
    h = 0.25
    while stab == false
        # Set up and down Kγ and KΩ values
        Kγu = (1+h)*Kγ
        Kγd = (1-h)*Kγ
        KΩu = (1+h)*KΩ
        KΩd = (1-h)*KΩ
        # Make corresponding microbes
        mγu = make_MicrobeP(MC,γm,ρ,Kγu,Pb,d,ϕH,KΩ,fd,R,Reacs,η,kcs,KSs,krs,n,ϕP)
        mγd = make_MicrobeP(MC,γm,ρ,Kγd,Pb,d,ϕH,KΩ,fd,R,Reacs,η,kcs,KSs,krs,n,ϕP)
        mΩu = make_MicrobeP(MC,γm,ρ,Kγ,Pb,d,ϕH,KΩu,fd,R,Reacs,η,kcs,KSs,krs,n,ϕP)
        mΩd = make_MicrobeP(MC,γm,ρ,Kγ,Pb,d,ϕH,KΩd,fd,R,Reacs,η,kcs,KSs,krs,n,ϕP)
        # Find corresponding χ2
        χγu = chi_square(mγu,S,P,ϕR0,Tmax,tsave,indm,admA,adsdA)
        χγd = chi_square(mγd,S,P,ϕR0,Tmax,tsave,indm,admA,adsdA)
        χΩu = chi_square(mΩu,S,P,ϕR0,Tmax,tsave,indm,admA,adsdA)
        χΩd = chi_square(mΩd,S,P,ϕR0,Tmax,tsave,indm,admA,adsdA)
        # print out data
        println("--------------------------------------------------------------")
        println("KΩ = $(KΩ), Kγ = $(Kγ), h = $(h)")
        println("Old chi squared = $(χ2)")
        println("Gamma up = $(χγu)")
        println("Gamma down = $(χγd)")
        println("Omega up = $(χΩu)")
        println("Omega down = $(χΩd)")
        # find smallest chi squared
        if χγu < χ2 && χγu < χγd && χγu < χΩu && χγu < χΩd
            println("Gamma up smallest")
            # Update chi squared to new smallest
            χ2 = χγu
            # Update relevant parameter
            Kγ = Kγu
        elseif χγd < χ2 && χγd < χΩu && χγd < χΩd
            println("Gamma down smallest")
            χ2 = χγd
            Kγ = Kγd
        elseif χΩu < χ2 && χΩu < χΩd
            println("Omega up smallest")
            χ2 = χΩu
            KΩ = KΩu
        elseif χΩd < χ2
            println("Omega down smallest")
            χ2 = χΩd
            KΩ = KΩd
        elseif h < 1e-5
            stab = true
        else
            println("Old value smallest")
            h /= 2.0
        end
    end
    # C and T now must be saved
    C, T = eval_a(m,S,P,ϕR0,Tmax,tsave)
    return(Kγ,KΩ,C,T)
end

# Need to work out how to remove the correct points
function atp_read()
    # Read in csv file
    dataf = DataFrame!(CSV.File("Data/dataset_22C_3d_atp_2.csv"))
    # Count the number of unique ID's
    ids = unique(dataf.ID)
    ni = length(ids)
    # count the number of trait names
    tns = unique(dataf.trait_name)
    nt = length(tns)
    # Need to make a data structure to contain all the data
    data = fill(Float64[],ni,nt,4,2)
    # Structure to store genus names
    genus = Array{String,1}(undef,ni)
    # Store genus names
    for i = 1:ni
        # Find all data of a certain type for a particular ID
        It = (dataf.ID .== ids[i])
        # Find genus name to use as a title
        genus[i] = dataf.bacterial_genus[It][1]
    end
    # Make a bool of valid points
    val = fill(true,length(dataf.ID))
    # Locate anamolous points
    an1 = (dataf.bacterial_genus .== "Flavobacterium") .& (dataf.replicate .== 2) .&
    (dataf.minute .>= 3012) .& (dataf.trait_name .== tns[3])
    an2 = (dataf.bacterial_genus .== "Pseudomonas(2)") .& (dataf.replicate .== 1) .&
    (dataf.minute .>= 3012) .& (dataf.trait_name .== tns[3])
    an3 = (dataf.bacterial_genus .== "Serratia") .& (dataf.replicate .== 1) .&
    (dataf.trait_name .== tns[3])
    # Set these points as invalid
    val[an1] .= false
    val[an2] .= false
    val[an3] .= false
    # Loop over all three conditions
    for m = 1:3
        for i = 1:ni
            # Find all cell count data with species ID
            I = (dataf.ID .== ids[i]) .& (dataf.trait_name .== tns[m])
            # Count number of replicates
            rps = unique(dataf.replicate[I])
            nr = length(rps)
            # Find unique times and order them
            ts = sort(unique(dataf.minute[I]))
            # Only remove points if they are ATP data
            if m == 3
                # Now loop over times
                for j = 1:length(ts)
                    # Group samples by time
                    IT = I .& (dataf.minute .== ts[j])
                    # Then perform hampel filter on sample
                    kp = hampel_filter(dataf.trait_value[IT],3)
                    # Loop over this output
                    for k = 1:length(kp)
                        # Overwrite the replicate number with 0
                        if kp[k] == 0
                            # Find indices of true values in IT
                            inds = findall(x->x==true,IT)
                            # The k'th index should then be marked false
                            val[inds[k]] = false
                        end
                    end
                end
            end
            # Loop over replicates
            for j = 1:nr
                # Only save valid cases with particular replicate number
                I2 = I .& (dataf.replicate .== rps[j]) .& val
                # Save times
                data[i,m,j,1] = dataf.minute[I2]
                # and save corresponding values
                if m == 3
                    data[i,m,j,2] = 1e-9.*dataf.trait_value[I2]
                else
                    data[i,m,j,2] = dataf.trait_value[I2]
                end
            end
        end
    end
    # Now want to find peak ATP values
    pt = zeros(ni)
    av = zeros(ni)
    for i = 1:ni
        # find all time points represented in ATP data
        ts = unique([data[i,3,1,1];data[i,3,2,1];data[i,3,3,1];data[i,3,4,1]])
        # loop over these time points
        for j = 1:length(ts)
            n = 0
            aT = 0.0
            # Loop over replicates
            for k = 1:4
                # Check if this time point is in the data
                if ts[j] ∈ data[i,3,k,1]
                    # increment count of points
                    n += 1
                    # Find index of point
                    ind = findfirst(x->x==ts[j],data[i,3,k,1])
                    # Add relevant point to the total
                    aT += (data[i,3,k,2])[ind]
                end
            end
            # Update values if average height is higher
            if aT/n > av[i]
                av[i] = aT/n
                if j > 1
                    pt[i] = ts[j-1]
                else
                    # Allow for negative start time
                    pt[i] = ts[1] - (ts[2]-ts[1])
                end
            end
        end
    end
    # Now find means and standard deviations of plots so that best fits can be found
    for i = 1#ni-1#1:ni
        # find all time points represented in ATP data
        ts = unique([data[i,3,1,1];data[i,3,2,1];data[i,3,3,1];data[i,3,4,1]])
        mA = zeros(length(ts))
        sdA = zeros(length(ts))
        # loop over time points
        for j = 1:length(ts)
            dataA = []
            n = 0
            for k = 1:4
                # Find if point it represented in the data
                if ts[j] ∈ data[i,3,k,1]
                    # If so increment counter
                    n += 1
                    # find relevant data point
                    ind = findfirst(x->x==ts[j],data[i,3,k,1])
                    # and add to the vector
                    dataA = cat(dataA,(data[i,3,k,2])[ind],dims=1)
                end
            end
            # Calculate mean
            mA[j] = sum(dataA)/n
            if n != 1
                # Find standard deviation in case where it can be calculated
                sdA[j] = stdm(dataA,mA[j])
            else
                # In case when mean calculated from one point, fix sd to high but not absurd value
                sdA[j] = 0.5*mA[j]
            end
        end
        # Give data to find best fit
        Kγ, KΩ, C, T = prot_fit(ts,pt[i],mA,sdA)
        # Plotting needs to change if there's an offset
        plot(T.+(pt[i]*60.0),C[:,1])
        plot!(ts*60.0,mA*6.02e23/1e9,yerror=sdA*6.02e23/1e9)
        savefig("Output/test.png")
        plot(T,C[:,1])
        savefig("Output/test1.png")
        plot(ts*60.0,mA*6.02e23/1e9)
        savefig("Output/test2.png")
        return(nothing)
    end
    # Call PyPlot
    pyplot()
    # Set a color-blind friendly palette
    theme(:wong2,dpi=150)
    # Plot all graphs
    for k = 1:3
        for i = 1:ni
            # Add appropriate times and labels
            plot(xlabel="Time (hours)",ylabel=tns[3],title=genus[i])
            # Loop over repeats to plot replicates
            for j = 1:4
                scatter!(data[i,k,j,1]/60.0,data[i,k,j,2],label="")
            end
            # Then save plot
            if k == 1
                savefig("Output/ATPDataPlots/Cells/$(genus[i]).png")
            elseif k == 2
                savefig("Output/ATPDataPlots/Biomass/$(genus[i]).png")
            else
                vline!([pt[i]/60.0],label="")
                savefig("Output/ATPDataPlots/ATP/$(genus[i]).png")
            end
        end
    end
    return(nothing)
end

@time atp_read()
