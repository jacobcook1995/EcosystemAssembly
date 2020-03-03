module Assembly

# Include scripts that contain other functions
include("MyPlots.jl") # This means I must always have a version of MyPlots.jl available
include("parameters/marsParameters.jl") # Parameters for base case from Marsland et al.
include("parameters/inhibParameters.jl") # Parameters for Inhibition case
include("simulate/mars.jl") # Include marsland model simulation code

# export global constants
export Rgas, ΔGATP

# Decleration of useful constants
global Rgas = 8.31446261815324 # gas constant in J.K^-1.mol^-1
global ΔGATP = 75000.0 # Need to find a reference for this at some point

end # module
