# Assembly Package
Made a package to contain code used for modelling microbial ecosystem assembly with the inclusion of end product inhibition.
The functions should be contained in a module in Assembly.jl.
I keep test scripts in test.jl.
Old functions that are no longer usable but worth saving are stored (temporarily) in grave.jl.
The scripts that I used to plot the graphs used in the Syntrophy SI document are saved in asmbplots.jl.
I have made a script called MyPlots.jl to store functions to assist plotting.
The simulate folder contains scripts that run simulations and return results.
I have a script called mars.jl that implements the model from Marsland et al.
In the same folder I also have a script called inhib.jl that implements an extended model with thermodynamic end product inhibition included.
The folder named parameters contains scripts that construct custom made Parameters Types.
The constructor for the parameter type used in the Marsland model is found in marsParameters.jl in the parameters folder.
The constructor for the parameter type used in our extended model is found in inhibParameters.jl in the parameters folder.
Have also created a folder called analysis where scripts analysis the output can be kept.
In this folder the script compare.jl attempts to compare output from the base and extended model for as similar as possible parameters.
I will have to decide what do when I have a lot of complicated functions at a later date.
At present all images are being plotted into an Output folder, this will have to be changed eventually.
All output data is stored in Data this can probably be preserved.\
