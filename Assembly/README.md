# Assembly Package
Made a package to contain code used for modelling microbial ecosystem assembly with the inclusion of end product inhibition.
The functions should be contained in a module in Assembly.jl.
I keep test scripts in test.jl.
Old functions that are no longer usable but worth saving are stored (temporarily) in grave.jl.
The scripts that I used to plot the graphs used in the Syntrophy SI document are saved in asmbplots.jl.
I have made a script called MyPlots.jl to store functions to assist plotting.
The folder named parameters contains scripts that construct custom made Parameters Types.
I have a script called mars.jl that implements the model from Marsland et al.
The constructor for the parameter type used in the Marsland model is found in marsParameters.jl in the parameters folder.
I will have to decide what do when I have a lot of complicated functions at a later date.
At present all images are being plotted into an Output folder, this will have to be changed eventually.
All output data is stored in Data this can probably be preserved.\
