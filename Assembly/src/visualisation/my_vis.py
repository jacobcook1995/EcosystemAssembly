# Edited most of these functions
import networkx as nx
import numpy as np
from os import listdir
import re
import importlib
import jonny_code
importlib.reload(jonny_code) # force reload, in case jonny_code.py is changed

# This is basically a massive function written by Jonny to wrangle the data
# big function to convert .csv into a networkx Graph object
# along with all the bells and whistles such as node/edge colors into attributes
# note that the functions in subsequent cells may break if attributes are removed
def wrangle(path, normalise_flux=False, width_min=1, width_max=5, alpha_min=.4, alpha_max=.6):
    indices = []
    with open(path) as f:
        for line in f.readlines():
            indices.append([])
            for index in line.split(','):
                indices[-1].append(index.strip())

    # # first set of targets should equal second set of sources
    # assert indices[1] == indices[2]
    # # interaction strengths are equal too
    # assert indices[5] == indices[6]

    # Find last forward slash in path so that the name can be extracted properly
    slh = path.rfind('/')
    # create nx graph
    G = nx.DiGraph()
    G.graph['name'] = path[slh+1:-4] # remove '.csv' and directories

    # local function to add edge to graph, sums flux if already in graph
    def add_flux(ij, flux):
        if G.has_edge(*ij):
            G.edges[ij]['weight'] += flux
        else:
            G.add_edge(*ij, weight=flux)

    # i is metabolite, j is species
    for i,j,w in zip(indices[0], indices[1], indices[4]):
        ij = ('m'+i, j)
        add_flux(ij, float(w))

    # i is species, j is metabolite
    for i,j,w in zip(indices[1], indices[2], indices[4]):
        ij = (i, 'm'+j)
        add_flux(ij, float(w))

    nx.set_node_attributes(G, {i: 1 if i in indices[1] else 0 for i in G.nodes}, 'x')
    nx.set_node_attributes(G, {i: int(i[1:]) if 'm' in i else int(i) for i in G.nodes}, 'y')

    # set node and edge colours (from matplotlib 'tab10')
    col_species = [.97,.50,.27]
    col_metabol = [.12,.47,.71]
    nx.set_node_attributes(G, {i:col_species if G.nodes[i]['x']==0 else col_metabol for i in G.nodes}, 'color')

    # set edge weight attributes according to flux on that interaction
    flux_logmax = np.log(max(nx.get_edge_attributes(G, 'weight').values()))
    flux_logmin = np.log(min(nx.get_edge_attributes(G, 'weight').values()))
    cols_edge = {}
    widths_edge = {}
    for ij in G.edges:
        flux_norm = (np.log(G.edges[ij]['weight'])-flux_logmin) / (flux_logmax-flux_logmin)
        if normalise_flux:
            G.edges[ij]['weight'] = flux_norm

        # set edge 'width' attribute for drawing later
        widths_edge[ij] = width_min + (width_max-width_min)*flux_norm

        # set edge 'color' attribute for drawing later
        alpha = alpha_min + (alpha_max-alpha_min)*flux_norm
        if G.nodes[ij[0]]['x']==1:
            cols_edge[ij] = col_metabol + [alpha]
        else:
            cols_edge[ij] = col_species + [alpha]

    nx.set_edge_attributes(G, widths_edge, 'width')
    nx.set_edge_attributes(G, cols_edge, 'color')

    G = nx.relabel.convert_node_labels_to_integers(G)
    return G

# Function to create a layered (bipartite) graph
def layered(path):
    G = wrangle(path, normalise_flux=True, alpha_min=.3, alpha_max=.4)
    x_constraint = { i:G.nodes[i]['x']   for i in G.nodes }

    ### NOTE: the following lines can be used/uncommented to additionally fix the y axis ###
    my = 4
    ### Just metabolites constrained
    y_constraint = { i:-G.nodes[i]['y']/2 for i in G.nodes if x_constraint[i]==0 }
    ### Just strains constrained
    # y_constraint = { i:G.nodes[i]['y']/2 for i in G.nodes if x_constraint[i]==1 }
    ### Metabolite and strains constrained
    # y_constraint = { i:G.nodes[i]['y']/2 for i in G.nodes }
    ### No constraint
    # y_constraint = None

    # weight_threshold can be used to cut weak links, but may break code if it disconnects the graph
    # weights run between 0.0 and 1.0, 0.0 causes no change 1.0 maximum change
    # wt = 0.3
    wt = 0.0
    jonny_code.stress(G, 'Output', y_constraint=y_constraint, x_constraint=x_constraint, weight_threshold=wt)

# This function seems to bundle the hierarchy
def bundle(path):
    G = wrangle(path, normalise_flux=False, alpha_min=.1, alpha_max=.7)
    # Change draw if I want to see the heirachy
    hierarchy = jonny_code.cluster(nx.Graph(G), draw=False)
    jonny_code.bundle(G, hierarchy, 'Output', merge_dist=.1, beta=.9, noderadius=.02)


# Call functions once
# saves layered drawing to 'figures/NetworkR=3rpt=37_stress.png'
# layered('Data/nets/R=3rpt=49.csv')
# Call functions multiple times
for path in listdir('Data/nets'):
    if re.match(r'.*\.csv$', path):
        layered(f'Data/nets/{path}')
        # bundle(f'Data/nets/{path}')
