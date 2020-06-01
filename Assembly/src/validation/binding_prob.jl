# A script to estimate the binding probability that we should use in our model
using Assembly
using LaTeXStrings
using DifferentialEquations
using Plots
import PyPlot

# function to calculate effective transciption rate
function wx(w::Float64,θx::Float64,a::Float64)
    ω = w*a/(θx+a)
    return(ω)
end

# function to calculate effective transciption rate with self inhibition
function wx(w::Float64,θx::Float64,a::Float64,Kq::Float64,hq::Int64,q::Float64)
    ω = (w*a/(θx+a))*(1/(1 + (q/Kq)^hq))
    return(ω)
end

function bind_dynamics!(dx::Array{Float64,1},x::Array{Float64,1},ps::ProtParameters,kb::Float64,
                        ku::Float64,dm::Float64,hq::Int64,Kq::Float64,w::Float64,θr::Float64,
                        θnr::Float64,wm::Float64,a::Float64,t::Float64)
    # First find the effective elongation rate
    γ = γs(a,ps)
    # Then find the growth rate
    λ = γ*(x[4]+x[5]+x[6])/ps.MC
    # First three equations are for mRNA
    dx[1] = wx(w,θr,a)
    dx[2] = wx(wm,θnr,a)
    dx[3] = wx(w,θnr,a,Kq,hq,x[9])
    for i = 1:3
        dx[i] += -kb*x[i]*x[7] + ku*x[i+3] + (γ*x[i+3])/ps.n[i] - dm*x[i] - λ*x[i]
    end
    # then 3 equations for the complexes
    for i = 4:6
        dx[i] = kb*x[i-3]*x[7] - ku*x[i] - (γ*x[i])/ps.n[i-3] - λ*x[i]
    end
    # Then equation for the free ribosomes
    for i = 7
        dx[i] = (γ*x[4])/ps.n[1] - λ*x[i]
        for j = 1:3
            dx[i] += (γ*x[j+3])/ps.n[j] - kb*x[j]*x[7] + ku*x[j+3]
        end
    end
    # Adding a two final equations for other protein fractions
    dx[8] = (γ*x[5])/ps.n[2] - λ*x[8]
    dx[9] = (γ*x[6])/ps.n[3] - λ*x[9]
    return(dx)
end

# Function to carry out the appropriate tests and then plot the results
function bind_prob()
    # mRNA-ribosome binding and unbinding are near the diffusion limit
    kb = (1.0/60)
    ku = (1.0/60)
    # This mRNA decay rate is taken from SI of Weiße et al. (Taniguchi Y, et al.)
    dm = (0.1/60)
    # Initialise parameter set for other parameters
    ps = initialise_prot(false)
    # Parameters for housekeeping protein auto inhibiton taken from Weiße et al (no other reference)
    hq = 4 # Steep autoihibition
    Kq = 150000.0/2.5 # Reducing this value
    # Set maximum mRNA synthesis rate, all considered the same for now
    # This is an approximation to Weiße et al where these parameters are optimised
    w = (900.0/60)
    wm = (4.0/60) # mRNA for enzymes is made slower apparently
    # Transcriptional thresholds scaled to match units of translational threshold
    θr = 400.0*(ps.Kγ/7.0)
    θnr = 4.0*(ps.Kγ/7.0)
    # Now sub the parameters in
    b_dyns!(dx,x,a,t) = bind_dynamics!(dx,x,ps,kb,ku,dm,hq,Kq,w,θr,θnr,wm,a,t)
    # Make simulation time span
    Tmax = 10000.0
    tspan = (0,Tmax)
    # Need to choose sensible initial conditions
    rn = 1.0*ones(3)
    # Initially no ribosomes are bound, this makes protein accounting simpler
    c = zeros(3)
    # sensible protein fraction
    xs = [0.275/ps.n[1],0.275/ps.n[2],0.45/ps.n[3]]*ps.MC
    # gather together as an initial condition
    x0 = [rn;c;xs]
    # Now need to choose the value of a to evaluate for
    a = 5e7
    # Then setup and solve the problem
    println("Simulation started.")
    prob = ODEProblem(b_dyns!,x0,tspan,a)
    sol = DifferentialEquations.solve(prob)
    pyplot(dpi=200)
    plot(sol.t,sol'[:,1:3])
    savefig("Output/RNA.png")
    plot(sol.t,sol'[:,4],label=L"c_r",xlabel="Time")
    plot!(sol.t,sol'[:,5],label=L"c_m")
    plot!(sol.t,sol'[:,6],label=L"c_q")
    plot!(sol.t,sol'[:,7],label="free")
    savefig("Output/bound.png")
    # Calculate fractions
    ϕR = (sol'[:,4].+sol'[:,5].+sol'[:,6].+sol'[:,7])*ps.n[1]/ps.MC
    ϕM = sol'[:,8]*ps.n[2]/ps.MC
    ϕQ = sol'[:,9]*ps.n[3]/ps.MC
    # Then plot them
    plot(sol.t,ϕR,label=L"\phi_R",xlabel="Time")
    plot!(sol.t,ϕM,label=L"\phi_M")
    plot!(sol.t,ϕQ,label=L"\phi_Q")
    plot!(sol.t,ϕR.+ϕM.+ϕQ,label=L"\phi_R+\phi_M+\phi_Q")
    savefig("Output/fractions.png")
    # Now find growth rate and plot
    γ = γs(a,ps)
    λ = γ*(sol'[:,4].+sol'[:,5].+sol'[:,6])/ps.MC
    plot(sol.t,λ,label="",xlabel="Time")
    savefig("Output/growth.png")
    return(nothing)
end

# Keep this for the time being but not sure that this model tells us much of interest
@time bind_prob()
