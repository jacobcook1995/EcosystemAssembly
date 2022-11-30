# Script to plot elements needed for figure 1
using TradeOff
using Plots
using LaTeXStrings
import PyPlot

# Function to plot the efficiency with changing ribosome fraction
function efficency_plot()
    # arbitrarily set the simulation type to 1
    sim_type = 1
    # Read in hard coded simulation parameters
    Np, Nt, M, d, μrange = sim_paras(sim_type)
    # Define the Reaction set
    Rs = [1, 3, 5, 7]
    # Set mratio as well
    mratio = 1e-2
    # Make a microbe so that the parameters are accessible
    mic = new_mic(M, Rs, d, μrange, mratio)
    # Preallocate a ribosome fractions to calculate this for
    ϕRs = collect(0.0:0.01:(1 - mic.ϕH))
    # Preallocate vector of maximum growth rates
    max_λ = zeros(length(ϕRs))
    # Loop over ribosome fractions
    for i in eachindex(ϕRs)
        # Find maximum growth rate (energy saturated)
        max_λ[i] = (mic.γm * mic.Pb / mic.n[1]) * ϕRs[i]
    end
    # Preallocate container to store the efficiencies
    effs = zeros(length(ϕRs))
    # Find efficiency for each step
    for i in eachindex(ϕRs)
        effs[i] = 1 / χs(ϕRs[i], mic)
    end
    # Setup plotting
    pyplot(dpi = 200, guidefontsize = 18, tickfontsize = 10)
    # Check if directory exists and if not make it
    if ~isdir("Output/Fig1")
        mkdir("Output/Fig1")
    end
    # Define latex labels
    phiR = L"\phi_R"
    chi = L"\chi"
    m1 = L"^{-1}"
    # Plot the two things
    plot(max_λ,
         effs,
         label = false,
         ylims = (0, Inf),
         xlabel = "Max growth rate (s$(m1))",
         linewidth = 2.5)
    plot!(ylabel = "Efficiency (aa ATP$(m1))")
    savefig("Output/Fig1/efficiency.png")
    return (nothing)
end

@time efficency_plot()
