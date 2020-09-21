# This script exists to extract network structures from the output of form_comm.jl
using Assembly
using JLD
using CSV
using DataFrames

# function to extract a bipartite metabolite-strain network
function bi_net()
    # Check that sufficent arguments have been provided
    if length(ARGS) < 2
        error("need to specify community and number of repeats")
    end
    # Preallocate the variables I want to extract from the input
    R = 0
    rpt = 0
    # Check that all arguments can be converted to integers
    try
        R = parse(Int64,ARGS[1])
        rpt = parse(Int64,ARGS[2])
    catch e
           error("both inputs must be integer")
    end
    # Check that simulation type is valid
    if R < 1
        error("each strain must have more than 1 reaction")
    end
    # Check that number of simulations is greater than 0
    if rpt < 1
        error("repeat number must be at least 1")
    end
    # First check that files exists
    pfile = "Data/Type$(R)/ParasType$(R)Run$(rpt).jld"
    if ~isfile(pfile)
        error("run $(rpt) is missing a parameter file")
    end
    ofile = "Data/Type$(R)/OutputType$(R)Run$(rpt).jld"
    if ~isfile(ofile)
        error("run $(rpt) is missing an output file")
    end
    # Then load in data
    ps = load(pfile,"ps")
    out = load(ofile,"out")
    # Make empty lists
    mp = Array{Int64,1}(undef,0) # Metabolite posistions
    ml = Array{Int64,1}(undef,0) # Metabolite links
    mf = Array{Float64,1}(undef,0) # Consumption fluxes
    sp = Array{Int64,1}(undef,0) # Strain posistions
    sl = Array{Int64,1}(undef,0) # Strain links
    sf = Array{Float64,1}(undef,0) # Production fluxes
    # Metabolite to strain link => consumption
    # Strain to metabolite link => production
    # Loop over all species
    for i = 1:ps.N
        # Loop over all reactions the species has
        for j = 1:ps.mics[i].R
            # Find indices of product + substrate
            indS = ps.reacs[ps.mics[i].Reacs[j]].Rct
            indP = ps.reacs[ps.mics[i].Reacs[j]].Prd
            # Find standard Gibbs free energy
            ΔG0 = ps.reacs[ps.mics[i].Reacs[j]].ΔG0
            # And η value
            η = ps.mics[i].η[j]
            # Find thermodynamic inhibition factor for this reaction
            θ1 = θ(out[ps.N+indS],out[ps.N+indP],ps.T,η,ΔG0)
            # Check that reaction actually runs
            if out[ps.N+indS] > 0.0 && θ1 < 0.99
                # Find enzyme as we need this for the flux
                E = Eα(out[2*ps.N+ps.M+i],ps.mics[i],j)
                # Find flux, same for both cases as we are looking at a one to one reaction
                q = qs(out[ps.N+indS],out[ps.N+indP],E::Float64,j,ps.mics[i],ps.T,ps.reacs[ps.mics[i].Reacs[j]])
                f = out[j]*q/NA
                # Add consumption link
                mp = cat(mp,indS,dims=1)
                ml = cat(ml,i,dims=1)
                mf = cat(mf,f,dims=1)
                # Add production link
                sp = cat(sp,i,dims=1)
                sl = cat(sl,indP,dims=1)
                sf = cat(sf,f,dims=1)
            end
        end
    end
    # Combine all 4 lists into one matrix to output
    A = zeros(Int64,4,length(mp))
    A[1,:] = mp
    A[2,:] = ml
    A[3,:] = sp
    A[4,:] = sl
    # Add 2 fluxes lists to also be output
    B = zeros(Float64,2,length(mp))
    B[1,:] = mf
    B[2,:] = sf
    # Then write out as a csv file
    CSV.write("Output/NetworkR=$(R)rpt=$(rpt).csv",DataFrame(A),writeheader=false)
    CSV.write("Output/NetworkR=$(R)rpt=$(rpt).csv",DataFrame(B),writeheader=false,append=true)
    return(nothing)
end

@time bi_net()
