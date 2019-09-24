module MyPlots

# This is a module to contain helpful scripts to aid my plotting

export annpos

# A function to return positions for labels
function annpos(datax::Array{Float64,1},datay::Array{Float64,1},δx=0.15::Float64,δy=0.125::Float64)
    # Need minimas and maximas
    xmax = maximum(datax)
    xmin = minimum(datax)
    # Calculate optimal x position for label
    posx = xmin - δx*(xmax-xmin)
    ymax = maximum(datay)
    ymin = minimum(datay)
    # Different formula for y as y axis extends a different length
    posy = ymax + δy*(ymax-ymin)
    return(posx,posy)
end

end
