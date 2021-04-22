import networkx as nx
import numpy as np
import scipy as sp
from collections import defaultdict
import matplotlib.pyplot as plt

def hashimoto_matrix(G):
    L = G.fresh_copy()
    # Create a graph specific edge function.
    get_edges = nx.generators.line._edge_func(G)
    for from_node in get_edges():
        # from_node is: (u,v) or (u,v,key)
        L.add_node(from_node)
        for to_node in get_edges(from_node[1]):
            # Modification here: make sure there aren't reciprocal edges.
            if to_node != (from_node[1], from_node[0]):
                L.add_edge(from_node, to_node)
    return L


alpha = 2.6
nodes = 20

while True:  
    s = np.array(nx.utils.powerlaw_sequence(nodes, alpha), dtype=np.int) #100 nodes, power-law exponent alpha - 1
    if np.sum(s)%2 == 0 and s[s == 0].shape[0] == 0:
        break

g = nx.configuration_model(list(s))
g = nx.Graph(g) # remove parallel edges
g.remove_edges_from(g.selfloop_edges())

gd = g.to_directed()

hashimoto_mat = hashimoto_matrix(gd)

L = nx.to_scipy_sparse_matrix(hashimoto_mat, nodelist=list(hashimoto_mat), dtype=float)
eigenvalue, eigenvector = sp.sparse.linalg.eigs(L.T, k=1, which='LR', maxiter=50, tol=0)

largest = np.abs(eigenvector.flatten().real)
norm = sp.sign(largest.sum()) * sp.linalg.norm(largest)
norml = largest/norm
# now find the nodes' centralities
ncenters = defaultdict(list)
for idx, n in enumerate(hashimoto_mat):
    # line graph nodes are edge tuples of original graph
    # want inedges by definition of centrality
    ncenters[n[1]].append(idx)


hcentrality = {}
for n in ncenters:
   hcentrality[n] = np.sum(norml[ncenters[n]])


eigenNodes = nx.eigenvector_centrality_numpy(g)
normen = np.linalg.norm(list(eigenNodes.itervalues()))

norm_en = {k: abs(v)*normen for k, v in eigenNodes.iteritems()}

pos = nx.spring_layout(g)

n_color = np.asarray([hcentrality[n] for n in g.nodes])
nx.draw(g, pos, nodelist=hcentrality.keys(), node_size=[v[1]*3000 for v in hcentrality.items()], node_color=n_color, cmap='viridis')

n_color = np.asarray([norm_en[n] for n in g.nodes])
nx.draw(g, pos, nodelist=norm_en.keys(), node_size=[v[1]*3000 for v in norm_en.items()], node_color=n_color, cmap='viridis')

