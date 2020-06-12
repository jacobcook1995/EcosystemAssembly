# Assembly Package
Made a package to contain code used for modelling microbial ecosystem assembly with the inclusion of end product inhibition.
The structure of the module is provided in Assembly.jl.
I have made a script called MyPlots.jl to store functions to assist plotting.\
The simulate folder contains scripts that run simulations and return results.
In this folder the script called inhib.jl that implements an extended model with thermodynamic end product inhibition included.
The script analytic.jl contains analytic functions to assist the simulation scripts.
The script proteome.jl runs a more complex model including proteome fractions for the single species (two metabolite) case.\
The folder named parameters contains scripts that construct custom made Parameters Types.
The constructor for the parameter type used in our extended model is found in inhibParameters.jl in the parameters folder.
The file make_inhib.jl provides functions to make various parameter sets for the model with inhibition.
The constructor for the parameter type used in our proteome fraction model is found in protParameters.jl in the parameters folder.
The file make_prot.jl provides functions to make various parameter sets for the proteome model.\
Have also created a folder called analysis where scripts that analyse the output can be kept.
In this folder the script compare.jl compares the output of the inhibition model for a number of different parameter choices.
The file atpdata.jl reads in a data file and plots the data as a scatter plot.
The file repeat.jl simulates an ecosystem with repeated invasions over a long time period.
The output from this can then be plotted using plot_repeat.jl.
Various tests of the single population proteome model are performed by singpop_prot.jl.\
At present all images are being plotted into an Output folder.
All input data is stored in Data.
Parameter sets are saved to the Paras folder.\
