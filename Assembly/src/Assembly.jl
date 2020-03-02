module Assembly

# Include scripts that contain other functions
include("MyPlots.jl") # This means I must always have a version of MyPlots.jl available
include("parameters/marsParameters.jl") # Parameters for base case from Marsland et al.
include("parameters/inhibParameters.jl") # Parameters for Inhibition case

end # module
