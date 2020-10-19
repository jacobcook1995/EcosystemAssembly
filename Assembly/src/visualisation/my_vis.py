# To write this I've copied Jonny's notebook, I need to work out what's going on so
# that I can figure out how to change things

# Think these are all about importing species python functions
import networkx as nx
import numpy as np
from os import listdir
import re

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

    # first set of targets should equal second set of sources
    assert indices[1] == indices[2]
    # interaction strengths are equal too
    assert indices[4] == indices[5]

    # create nx graph
    G = nx.DiGraph()
    G.graph['name'] = path[:-4] # remove '.csv'

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
    for i,j,w in zip(indices[2], indices[3], indices[5]):
        ij = (i, 'm'+j)
        add_flux(ij, float(w))

    nx.set_node_attributes(G, {i: 0 if i in indices[1] else 1 for i in G.nodes}, 'y')
    nx.set_node_attributes(G, {i: int(i[1:]) if 'm' in i else int(i) for i in G.nodes}, 'x')

    # set node and edge colours (from matplotlib 'tab10')
    col_species = [.97,.50,.27]
    col_metabol = [.12,.47,.71]
    nx.set_node_attributes(G, {i:col_species if G.nodes[i]['y']==1 else col_metabol for i in G.nodes}, 'color')

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
        if G.nodes[ij[0]]['y']==0:
            cols_edge[ij] = col_metabol + [alpha]
        else:
            cols_edge[ij] = col_species + [alpha]

    nx.set_edge_attributes(G, widths_edge, 'width')
    nx.set_edge_attributes(G, cols_edge, 'color')

    G = nx.relabel.convert_node_labels_to_integers(G)
    return G
