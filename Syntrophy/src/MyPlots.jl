# This is a script that contains a number of functions to help with plotting
# Modules can include this script to gain access to the exported functions
using LsqFit
# THIS IS MY BEST GUESS FOR WHERE TO PUT THIS BUT MIGHT NEED TO BE MOVED IN FUTURE

# Export functions that are useful externally
export annpos corrparr

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

# Function that takes sets of data and does some fitting
# This returns the fit, the correlation coefficient, and the confidence interval
function corrparr(xdata::Array{Float64,1},ydata::Array{Float64,1},p0::Array{Float64,1},xran::AbstractRange)
    # p0 is the initial parameter guess
    # xran is the range to calculate over

    # First define model
    @. model(x,p) = p[1] + p[2]*x
    # Fit data to this model
    fit = curve_fit(model,xdata,ydata,p0)
    yint = coef(fit)[1]
    slop = coef(fit)[2]
    # Then find confidence interval
    con = confidence_interval(fit, 0.05)
    vyint = con[1]
    vslop = con[2]
    # Calculate intervals for the range given
    intlow = abs.(model(xran,[yint,slop]) .- model(xran,[vyint[2],vslop[2]]))
    intup = abs.(model(xran,[yint,slop]) .- model(xran,[vyint[1],vslop[1]]))
    # Now calculate Pearson correlation coefficient
    xbar = sum(xdata)/length(xdata)
    ybar = sum(ydata)/length(ydata)
    a = 0
    b = 0
    c = 0
    for i = 1:length(xdata)
        a += (xdata[i] - xbar)*(ydata[i] - ybar)
        b += (xdata[i] - xbar)^2
        c += (ydata[i] - ybar)^2
    end
    r = a/sqrt(b*c)
    return(yint,slop,intlow,intup,r)
end
