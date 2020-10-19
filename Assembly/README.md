# Assembly Package
Made a package to contain code used for modelling microbial ecosystem assembly with the inclusion of end product inhibition.
The structure of the module is provided in Assembly.jl.
I have made a script called MyPlots.jl to store functions to assist plotting.\
The simulate folder contains scripts that run simulations and return results.
In this folder the script called full.jl that implements the full model with thermodynamic end product inhibition and proteome constraints included.
The script analytic.jl contains analytic functions to assist the simulation scripts.
The script proteome.jl runs a simple version of the proteome constraint model for the single species (two metabolite) case.\
The folder named parameters contains scripts that construct custom made Parameters Types.
The constructor for the parameter type used in our proteome fraction model is found in protParameters.jl in the parameters folder.
The file make_prot.jl provides functions to make various parameter sets for the proteome model.
The constructor for the parameter type used in our full model is found in sullParameters.jl in the parameters folder.
The file make_full.jl provides functions to make various parameter sets for the full model.\
Have also created a folder called analysis where scripts that analyse the output can be kept.
The file atpdata.jl reads in a data file and plots the data as a scatter plot, fits to this data are then found using a gradient descent method.
The file repeat.jl simulates an ecosystem with repeated invasions over a long time period, it's currently redundant (i.e. doesn't run) and is just being saved as an example for future work.
Simulations of the full model are run using form_comm.jl, the output is then plotted using plot_comm.jl.
This data can be used to generate networks through the script net.jl.
The script testing.jl allows me to check why particular simulations crashed or otherwise failed.
Various tests of the single population proteome model are performed by singpop_prot.jl.\
The folder theory is for general theoretical results that are needed for the paper but not the simulation.
The script 2step.jl investigates the impact on the reaction thermodynamics of making the overall reaction multistep.\
The folder visualisation contains code for visualising microbial networks, the code in it was orginally written by Dr Jonathan X Zheng.
This code was adapted from his work in doi.org/10.1109/TVCG.2019.2944619 and doi.org/10.1109/TVCG.2018.2859997.\
At present all images are being plotted into an Output folder.
All input data is stored in Data.
Parameter sets are saved to the Paras folder.\
