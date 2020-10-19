# Thes function were all written by Jonny and I have (as of yet) made no alterations to them

import networkx as nx
import matplotlib.pyplot as plt
import scipy
import scipy.cluster
import numpy as np
import math
import cairosvg

###############
# stress layout

def stress(G, output_folder, x_constraint=None, y_constraint=None, weight_threshold=0):
    # remove weak edges
    G2 = G.copy()
    for ij in G.edges:
        if G.edges[ij]['weight'] < weight_threshold:
            G2.remove_edge(*ij)

    # compute layout
    X = _sgd(G2, y_constraint=y_constraint, x_constraint=x_constraint)

    # draw with colours
    cols_node = list(nx.get_node_attributes(G, 'color').values())
    cols_edge = list(nx.get_edge_attributes(G, 'color').values())
    widths_edge = list(nx.get_edge_attributes(G, 'width').values())
    labels = nx.get_node_attributes(G, 'x')
    nx.draw(G, pos=X, node_color=cols_node, edge_color=cols_edge, width=widths_edge, labels=labels, arrows=False)
    plt.axis('equal')
    plt.savefig(f'{output_folder}/{G.graph["name"]}_stress.png')

def _sgd(G: nx.Graph, x_constraint=None, y_constraint=None, t_max=30, eps=.1, t_min=-10, mu_max=1.1):
    # make undirected
    S = nx.to_scipy_sparse_matrix(G)
    d = scipy.sparse.csgraph.shortest_path(S, directed=False, unweighted=True)

    # initialise terms
    n = len(G)
    m = (n*(n-1))//2
    terms = np.empty((m, 4))
    ij = 0
    for i in range(n):
        for j in range(i):
            terms[ij] = [i,j,d[i,j],1/d[i,j]**2]
            ij += 1

    # set step sizes
    eta_min = eps / np.max(terms, axis=0)[3]
    eta_max = 1 / np.min(terms, axis=0)[3]
    decay = np.log(eta_min / eta_max) / (t_max - 1)
    schedule = eta_max * np.exp(decay * np.arange(t_min,t_max)) # t_min < 0 means extra init iterations

    # initialise positions
    X = np.random.random((n,2))

    # multiplier array to fix node positions
    pins = np.ones((n,2))

    if x_constraint is not None:
        for i,x in x_constraint.items():
            pins[i,0] = 0
            X[i,0] = x
    if y_constraint is not None:
        for i,y in y_constraint.items():
            pins[i,1] = 0
            X[i,1] = y

    # do optimisation
    for eta in schedule:
        np.random.shuffle(terms)
        for term in terms:
            i,j,d,w = int(term[0]), int(term[1]), term[2], term[3]
            mu = min(eta*w, mu_max)

            X_ij = X[i] - X[j]
            mag = np.linalg.norm(X_ij)
            r = (mag-d)/2 * X_ij/mag

            # the movements will be multiplied by 1|0 such to pin nodes
            X[i] -= mu * r * pins[i]
            X[j] += mu * r * pins[j]

    return X


###############################################
# clustering as preprocessing step for bundling

def cluster(G, steps=5, draw=False):
    embedding = _embed_snapshot_markov(G, steps, False)
    hierarchy = _cluster_walktrap(embedding, opt_ord=True)
    if draw:
        _draw_dendrogram(G, hierarchy)
    return scipy.cluster.hierarchy.to_tree(hierarchy)

def _embed_snapshot_markov(G, steps, degrees):
    '''embedding of pons & lapaty, including extra weighting using degrees'''
    # get markov matrix P
    n = len(G)
    A = nx.to_numpy_array(G)
    D = np.identity(n) * np.sum(A, axis=1)

    P = np.linalg.inv(D) @ A

    # steps
    P = np.linalg.matrix_power(P, steps)
    if degrees:
        P = np.linalg.inv(sqrtm(D)) @ P
    return P

def _cluster_walktrap(embedding, method='ward', opt_ord=False):
    dist = scipy.spatial.distance.pdist(embedding, metric='euclidean')
    return scipy.cluster.hierarchy.linkage(dist, method=method, optimal_ordering=opt_ord)

def _draw_dendrogram(G, hierarchy, figsize=(10,5), nodelabels=None):
    plt.figure(figsize=figsize)
    plt.rcParams["font.weight"] = "bold"
    dend = scipy.cluster.hierarchy.dendrogram(hierarchy, leaf_rotation=90, leaf_font_size=8)

    ticks, labels = plt.xticks()
    for label in labels:
        idx = int(label.get_text())
        label.set_color(G.nodes[idx]['color'])

    if nodelabels is not None:
        labels = [nodelabels[int(l.get_text())] for l in labels]
    plt.xticks(ticks, labels)
    plt.show()

######################################
# bundling with preprocessed hierarchy

def bundle(G, clustroot, output_folder, merge_dist=.1, tutte_expo=3, skip_lca=True, skip_intra=False, beta=.95, linkwidth=.005, linkopacity=.25, noderadius=.012):
    root = _prune(clustroot, merge_dist, tutte_expo)
    _route(root, G, skip_lca, skip_intra, plot_tree=False)
    svg = _draw_svg(G, beta=beta, linkwidth=linkwidth, linkopacity=linkopacity, noderadius=noderadius)
    cairosvg.svg2png(bytestring=svg, write_to=f'{output_folder}/{G.graph["name"]}_bundled.png', dpi=400)

class TreeNode:
    '''needed because ClusterNode for dendrogram is binary'''
    def __init__(self, idx, dist):
        self.idx = idx
        self.dist = dist
        self.children = []

        self.parent = None
        self.length = None
        self.angle = None

def _prune(clustroot: scipy.cluster.hierarchy.ClusterNode, mingap, length_exp=1) -> TreeNode:
    # first build tree with merged branches
    def build(node: scipy.cluster.hierarchy.ClusterNode):
        if node is None:
            return None

        here = TreeNode(node.id, node.dist/clustroot.dist) # normalise
        left = build(node.left)
        right = build(node.right)

        # if leaf, just return
        if left is None:
            assert right is None # should be satisfied from dendrogram
            return here

        # merge child into self if not leaf and gap is smaller than min
        if len(left.children) != 0 and here.dist-left.dist < mingap:
            here.children.extend(left.children)
        else:
            here.children.append(left)
        # and same for right
        if len(right.children) != 0 and here.dist-right.dist < mingap:
            here.children.extend(right.children)
        else:
            here.children.append(right)
        return here

    root = build(clustroot)

    # count number of leaves and set parent pointers
    nleaves = 0
    nbranches = 0
    def count(node: TreeNode):
        if node is None:
            return
        if len(node.children) == 0:
            nonlocal nleaves
            nleaves += 1
        else:
            nonlocal nbranches
            nbranches += 1
            for child in node.children:
                child.parent = node
                count(child)

    count(root)

    # place each node at the midangle of its children
    n = nleaves + nbranches
    leafcounter = 0
    branchcounter = 0
    def place(node: TreeNode):
        assert node is not None
        if len(node.children) == 0:
            nonlocal leafcounter
            node.angle = (leafcounter / nleaves) * 2 * math.pi
            node.length = 1
            leafcounter += 1
            return 1
        else:
            total_leaves = sum(place(child) for child in node.children)
            node.length = (1 - (total_leaves / nleaves)) ** length_exp # exponent makes layout more uniform
#             node.angle = (node.children[0].angle + node.children[-1].angle) / 2
            node.angle = sum(c.angle for c in node.children) / len(node.children)

            # reassign indices to make contiguous
            nonlocal branchcounter
            node.idx = nleaves + branchcounter
            branchcounter += 1
            return total_leaves

    place(root)
    return root

def _route(root: TreeNode, G: nx.Graph, removeLCA=True, remove_intra=False, plot_tree=False):
    # first convert to cartesian coordinates
    idx2node = {}
    node2pos = {}
    def cartesian(node):
        if node is None:
            return

        x = node.length * math.cos(node.angle)
        y = node.length * math.sin(node.angle)
        node2pos[node] = np.array([x, y])
        idx2node[node.idx] = node
        for child in node.children:
            cartesian(child)

    cartesian(root)
    idx2pos = { i:node2pos[x] for i,x in idx2node.items() }
    nx.set_node_attributes(G, idx2pos, 'pos')

    colors = nx.get_node_attributes(G, 'color')

    if plot_tree:
        def drawtree(node):
            if node is None:
                return
            x0,y0 = node2pos[node]
            if node.parent is not None:
                x1,y1 = node2pos[node.parent]
                ax.plot([x0,x1],[-y0,-y1], color='k')
            for child in node.children:
                drawtree(child)
        fig, ax = plt.subplots()
        drawtree(root)
        ax.axis('equal')
        plt.show()

    def LCA(src1, src2):
        curr1 = src1
        curr2 = src2
        while curr1 is not curr2:
            curr1 = curr1.parent if curr1 is not None else src2
            curr2 = curr2.parent if curr2 is not None else src1
        return curr1

    def getpath(src, tgt):
        path = []
        while src is not tgt:
            path.append(node2pos[src])
            src = src.parent
        return path

    # then route edges through hierarchy
    splines = {}
    for src,tgt in G.edges:
        lca = LCA(idx2node[src], idx2node[tgt])
        assert lca is not None

        src2lca = getpath(idx2node[src], lca)
        tgt2lca = getpath(idx2node[tgt], lca)

        if not removeLCA or len(src2lca)+len(tgt2lca) <= 2:
            if not remove_intra:
                splines[(src,tgt)] = np.array(src2lca + [node2pos[lca]] + list(reversed(tgt2lca)))
        else:
            splines[(src,tgt)] = np.array(src2lca + list(reversed(tgt2lca)))

    nx.set_edge_attributes(G, splines, 'spline')

def draw_bspline_cubic(path, color=None) -> str:
    """draws a cubic b-spline, with an open knot vector and repeated control points"""
    svg = []
    m = len(path)
    if m < 2:
        raise ValueError("path is less than 2 points long")
    elif m == 2:
        p0 = path[0]
        p1 = path[1]
        svg.append(f'<path d="M {p0[0]:.1f} {p0[1]:.1f} L {p1[0]:.1f} {p1[1]:.1f}"/>')
    else:
        p000 = path[0] # not strictly correct, but works
        p112 = 2/3*path[0] + 1/3*path[1]
        p122 = 1/3*path[0] + 2/3*path[1]
        svg.append(f'<path d="M {p000[0]:.1f} {p000[1]:.1f} C {p112[0]:.1f} {p112[1]:.1f} {p122[0]:.1f} {p122[1]:.1f}')

        for i in range(1, len(path)-1):
            p123 = path[i]
            p234 = path[i+1]
            p223 = 2/3*p123 + 1/3*p234
            p233 = 1/3*p123 + 2/3*p234
            p222 = .5*p122 + .5*p223

            svg.append(f' {p222[0]:.1f} {p222[1]:.1f} S {p233[0]:.1f} {p233[1]:.1f}')
            p122 = p233

        end = path[-1]
        svg.append(f' {end[0]:.1f} {end[1]:.1f}"')

        if color:
            svg.append(f' stroke="rgb({int(color[0]*255)},{int(color[1]*255)},{int(color[2]*255)})" stroke-opacity="{color[3]/2:.3f}"')
        svg.append(f'/>')

    return(''.join(svg))

def _draw_svg(G, beta=.75, width=750, border=50, linkwidth=.05, noderadius=.1, linkopacity=1, nodeopacity=1) -> str:
    scale = width-2*border
    svg = []
    svg.append(f'<svg width="{width:.0f}" height="{width:.0f}" xmlns="http://www.w3.org/2000/svg">')
    svg.append('<style type="text/css">')
#     svg.append(f'path{{stroke:black;stroke-width:{scale*linkwidth:.3f};stroke-opacity:{linkopacity:.3f};stroke-linecap:round;fill:transparent}}')
    svg.append(f'path{{stroke-width:{scale*linkwidth:.3f};stroke-linecap:round;fill:transparent}}')
    svg.append(f'circle{{r:{scale*noderadius:.3f};stroke-width:0;fill-opacity:{nodeopacity:.3f}}}')
    svg.append('</style>')

    # add white background
    svg.append('<rect width="100%" height="100%" fill="white"/>')

    # draw splines
    colors = nx.get_edge_attributes(G, 'color')
    splines = nx.get_edge_attributes(G, 'spline')
    for edge,spline in splines.items():
        # rescale
        for i in range(len(spline)):
            spline[i] = [border,border] + (spline[i]+[1,1])/2*scale
        # straighten with beta
        for i in range(1,len(spline)-1):
            spline[i] = beta*spline[i] + (1-beta) * (spline[0] + (i/(len(spline-1)))*(spline[-1]-spline[0]))

        svg.append(draw_bspline_cubic(spline, colors[edge]))

    colors = nx.get_node_attributes(G, 'color')
    positions = nx.get_node_attributes(G, 'pos')
    labels = nx.get_node_attributes(G, 'x')
    for node in colors.keys():
        pos = [border,border] + (positions[node]+[1,1])/2*scale
        color = f'rgb({int(colors[node][0]*255)},{int(colors[node][1]*255)},{int(colors[node][2]*255)})'
        svg.append(f'<circle cx="{pos[0]:.1f}" cy="{pos[1]:.1f}" fill="{color}"/>')
        svg.append(f'<text x="{pos[0]:.1f}" y="{pos[1]:.1f}" font-size="{width/50}" text-anchor="middle" alignment-baseline="central">{labels[node]}</text>')

    svg.append('</svg>')

    return '\n'.join(svg)
