(('micro_mg1_0.F90', 'micro_mg1_0', 'micro_mg_tend'), '0.3995832581522189'), (('wv_saturation.F90', 'wv_saturation', 'svp_water'), '0.0153547485783434'), (('wv_sat_methods.F90', 'wv_sat_methods', 'wv_sat_svp_water'), '0.0017542166910599058'), (('wv_sat_methods.F90', 'wv_sat_methods', 'goffgratch_svp_water_r4'), '0.00020041202130351862'), (('wv_sat_methods.F90', 'wv_sat_methods', 'wv_sat_qsat_water'), '7.951300041923211e-06'), (('wv_sat_methods.F90', 'wv_sat_methods', 'wv_sat_svp_to_qsat'), '9.084032329122571e-07')


nx.draw(weak, pos, node_color='#3778bf', edge_color=colors,
         width=1.25, edge_cmap=plt.cm.Greens, with_labels=False, node_size=100)



labels = {}    
for node in weak.nodes:
    if 'goffgratch_svp_water_r4' in node:
        labels[node] = 'GoffGratch_svp_water'
        break

sizes = []
colors = []
for node in weak.nodes:
    if 'goffgratch_svp_water_r4' in node:
        sizes.append(300)
	colors.append('#cf000f')
    else:
        sizes.append(160)
	colors.append('#3778bf')      

nx.draw_networkx_labels(weak,pos,labels,font_size=24,font_color='black')



eqRel = lambda u, v: (CESM.node[u]['attrs']['loc'][0] == CESM.node[v]['attrs']['loc'][0])

posq = nx.spring_layout(q)
posqnew = copy.deepcopy(posq)

labelsq = {}    
for node in q.nodes:
    name = CESM.node[iter(node).next()]['attrs']['loc'][0]
    print name
    if 'wv_sat' in name:
        labelsq[node] = name
        posqnew[node] = posqnew[node] + np.array([0., 0.08])



sizesq = []
colorsq = []
for node in q.nodes:
    name = CESM.node[iter(node).next()]['attrs']['loc'][0]
    if 'wv_sat' in name:
        sizesq.append(600)
	colorsq.append('#cf000f')
    else:
        sizesq.append(300)
	colorsq.append('#3778bf')

nx.draw(q, posq, node_color=colorsq, edge_color='#7bb274',
        width=3, with_labels=False, node_size=sizesq)

nx.draw_networkx_labels(q,posqnew,labelsq,font_size=34,font_color="black")


>>> gg4 = set() 
>>> for n in e:
...   if 'goffgratch_svp_water_r4' in n:
...     gg4.add(n)

gg4sg = nx.subgraph(CESM.graph, e)

innodes = set()
for n in f[:10]:
  for p in nx.shortest_path(gg4sg, target=n[0]):
      innodes.add(p)

insubgraph = nx.subgraph(gg4sg, innodes)

pos = nx.spring_layout(insubgraph)

central = set([i[0] for i in f[:10]])

sizes = []
colors = []
for node in insubgraph.nodes:
    if node in gg4:
        sizes.append(600)
        colors.append('#cf000f')
    elif node in central:
        sizes.append(800)
        colors.append('#fd6a02')
    else:
        sizes.append(100)
        colors.append('#3778bf')


pathset = set()
for n in f[:10]:
    for bad in gg4:
        if nx.has_path(insubgraph, bad, n[0]):
            for p in pairwise(nx.shortest_path(insubgraph, source=bad, target=n[0])):
                pathset.add(p)

edgecolor = []
edgewidth = []
for edge in insubgraph.edges:
    if edge in pathset:
        edgecolor.append('#60526C')
        edgewidth.append(3)
    else:
        edgecolor.append('#7bb274')
        edgewidth.append(0.5)

nx.draw(insubgraph, pos, node_color=colors, edge_color=edgecolor,
        width=edgewidth, with_labels=False, node_size=sizes, alpha=0.8)

>>> for n in f[:400]:
...   if nx.has_path(dyn3, 'omega_p__compute_and_apply_rhs', n[0]):
...     print n[0]
... 
elem__compute_and_apply_rhs
elem__euler_step

>>> innodes = set()
>>> for n in f[:100]:
...   for p in nx.shortest_path(dyn3, target=n[0]):
...       innodes.add(p)





import matplotlib.pyplot as plt
import seaborn as sns
from metagraph_utils import *
import networkx as nx
from subgraphs import *
def plot_communities(subgraph, comms, sourcest, thresh, numCenters, name, layout='spring_layout', centerColor='#45ADA8', badColor='#cf000f', pos=None, pdf=False):
    sources = set()
    for n in subgraph:
        if isinstance(sourcest, str):
            if sourcest == n:
                sources.add(n)
        elif isinstance(sourcest, list):
            for sstr in sourcest:
                if sstr == n:
                    sources.add(n)
    print 'number of sources: ', sources
    if not pos:
        if layout == 'spring_layout':
            pos = nx.spring_layout(subgraph)
        elif layout == 'spectral_layout':
            pos = nx.spectral_layout(subgraph)
        elif layout == 'random_layout':
            pos = nx.random_layout(subgraph)
    nodesizes = []
    nodecolors = []
    for n in subgraph.nodes:
        nodecolor = sns.color_palette("Paired").as_hex()[1]
        size = 25
        if n in sources:
            nodecolor = badColor
            size = 300
        nodecolors.append(nodecolor)
        nodesizes.append(size)
    plt.figure(figsize=(16, 9), dpi=200)
    nx.draw(subgraph, pos, node_color=nodecolors, edge_color=sns.color_palette("Paired").as_hex()[0],
        width=0.5, with_labels=False, node_size=nodesizes, alpha=0.85)
    if not pdf:
        plt.savefig(name + '_nocomms' + '.png', format='png')
    else:
        plt.savefig(name + '_nocomms' + '.pdf', format='pdf')
    basenode = '#000000'
    baseedge = '#000000'
    bigcomms = [c for c in comms if len(c['subgraph']) > 1]
    pal = sns.color_palette("Paired", 2*len(bigcomms)).as_hex()
    colormap = []
    sourceidx = None
    for idx, c in enumerate(bigcomms):
        tmp = {}
        tmp['nodes'] = set(c['subgraph'].nodes)
        tmp['edges'] = set(c['subgraph'].edges)
        tmp['ncolors'] = pal[2*idx + 1]
        tmp['edgecolors'] = pal[2*idx]
        colormap.append(tmp)
        if sources & tmp['nodes']:
            sourceidx = idx
    nodesizes = []
    nodecolors = []
    for n in subgraph.nodes:
        nodecolor = '#000000'
        size = 25
        for cn in colormap:
            if n in cn['nodes']:
                nodecolor = cn['ncolors']
                edgecolor = cn['edgecolors']
        if n in sources:
            nodecolor = badColor
            size = 300
        nodecolors.append(nodecolor)
        nodesizes.append(size)
    edgecolors = []
    for edg in subgraph.edges:
        edgecolor = '#000000'
        for cn in colormap:
            if edg in cn['edges']:
                edgecolor = cn['edgecolors']
        edgecolors.append(edgecolor)
    plt.figure(figsize=(16, 9), dpi=200)
    nx.draw(subgraph, pos, node_color=nodecolors, edge_color=edgecolors,
        width=0.5, with_labels=False, node_size=nodesizes, alpha=0.85)
    if not pdf:
        plt.savefig(name + '.png', format='png')
    else:
        plt.savefig(name + '.pdf', format='pdf')
    srcsubgraph = bigcomms[sourceidx]
    srcpos = nx.spring_layout(srcsubgraph['subgraph'])
    pathset = set()
    centers = set()
    center = 1.
    count = min(numCenters, len(srcsubgraph['centralities']))
    for i in xrange(count):#xrange(len(srcsubgraph['centralities'])):
        center = srcsubgraph['centralities'][i][1]
        if center < thresh:
            break           
        n = srcsubgraph['centralities'][i][0]
        centers.add(n)
        if n in sources:
            print 'found'
        for bad in sources:
            if bad in srcsubgraph['subgraph'] and n in srcsubgraph['subgraph']:
                if nx.has_path(srcsubgraph['subgraph'], bad, n):
                    for p in pairwise(nx.shortest_path(srcsubgraph['subgraph'], source=bad, target=n)):
                        pathset.add(p)
    nodesizes = []
    nodecolors = []
    edgecolors = []
    edgewidths = []
    edgelist = []
#    alphal = []
    for sn in srcsubgraph['subgraph'].nodes:
        nodecolor = colormap[sourceidx]['ncolors']
        size = 25
        if sn in sources and sn in centers:
            nodecolor = badColor
            size = 700            
        elif sn in sources:
            nodecolor = badColor
            size = 300
        elif sn in centers:
            nodecolor = centerColor
            size = 400
        nodecolors.append(nodecolor)
        nodesizes.append(size)
    for edg in srcsubgraph['subgraph'].edges:
        edgecolor = colormap[sourceidx]['edgecolors']
        width = 0.5
        if edg not in pathset:
#            edgecolor = '#60526C'
#            width = 3
            edgecolors.append(edgecolor)
            edgewidths.append(width)
            edgelist.append(edg)
    # for edge in pathset:
    #     edgecolors.append('#60526C')
    #     edgewidths.append(3)
    #     edgelist.append(edge)
    #     alphal.append(1.0)
    plt.figure(figsize=(16, 9), dpi=200)
    nx.draw_networkx_nodes(srcsubgraph['subgraph'], pos, node_color=nodecolors,
        node_size=nodesizes, alpha=0.85)   
    nx.draw_networkx_edges(srcsubgraph['subgraph'], pos, edgelist=edgelist, edge_color=edgecolors,
        width=edgewidths, alpha=0.85)
    #nx.draw(srcsubgraph['subgraph'], pos, node_color=nodecolors, edge_color=edgecolors,
    #    width=edgewidths, with_labels=False, node_size=nodesizes, alpha=0.85)
    if pathset:
        nx.draw_networkx_edges(srcsubgraph['subgraph'], pos, edgelist=pathset, alpha=1.0, width=3, edge_color='#60526C')
    plt.axis('off')
    if not pdf:
        plt.savefig(name + '_srcsubgraph' + '.png', format='png')
    else:
        plt.savefig(name + '_srcsubgraph' + '.pdf', format='pdf')
    newsg = None
    if not pathset:
        ysgn = set()
        notsgn = set()
        for bc in bigcomms:
            ssg = bc['subgraph']
            count = min(numCenters, len(bc['centralities']))
            for i in xrange(count):#xrange(len(bc['centralities'])):
                center = bc['centralities'][i][1]
                if center < thresh:
                    break
                n = bc['centralities'][i][0]
                print n
                for bad in sources:
                    tf = nx.has_path(subgraph, bad, n)
                    for nn in nx.shortest_path(subgraph, target=n):
                        if tf:
                            ysgn.add(nn)
                        else:
                            notsgn.add(nn)
        if ysgn:
            newsg = nx.subgraph(subgraph, ysgn)
            print 'connected to other'
        else:
            notsg = set(subgraph) - notsgn
            newsg = nx.subgraph(subgraph, notsg)
            print 'not connected', len(notsgn)
    else:
        newsgn = set()
        newsgc = set()
        for c in centers:
            for bad in sources:
                if bad in srcsubgraph['subgraph']:
                    if nx.has_path(srcsubgraph['subgraph'], bad, c):
                        newsgc.add(c)
                        newsgn.add(c)
        for nc in newsgc:
            for n in nx.shortest_path(subgraph, target=nc):
                newsgn.add(n)
        newsg = nx.subgraph(subgraph, newsgn)
        print 'connected', len(newsg)
    return pos, newsg

from math import log10
import collections
def drop_zeros(a_list):
    return [i for i in a_list if i>0]

def log_binning(counter_dict,bin_count=35):
    max_x = log10(max(counter_dict.keys()))
    max_y = log10(max(counter_dict.values()))
    max_base = max([max_x,max_y])
    min_x = log10(min(drop_zeros(counter_dict.keys())))
    bins = np.logspace(min_x,max_base,num=bin_count)
    # Based off of: http://stackoverflow.com/questions/6163334/binning-data-in-python-with-scipy-numpy
    bin_means_y = (np.histogram(counter_dict.keys(),bins,weights=counter_dict.values())[0] / np.histogram(counter_dict.keys(),bins)[0])
    bin_means_x = (np.histogram(counter_dict.keys(),bins,weights=counter_dict.keys())[0] / np.histogram(counter_dict.keys(),bins)[0])
    return bin_means_x,bin_means_y


degree_sequence = sorted([d for n, d in CESM.graph.degree()], reverse=True)  # degree sequence
# print "Degree sequence", degree_sequence
degreeCount = collections.Counter(degree_sequence)
deg, cnt = zip(*degreeCount.items())

cesm_x, cesm_y = log_binning(degreeCount, 100)
plt.figure(figsize=(10.666666667, 6), dpi=200)
plt.xscale('log')
plt.yscale('log')
plt.scatter(deg, cnt, c='#6C5B7B')
plt.xlim((0.9, 1200))
plt.ylim((0.9, 35000))
#plt.axis('tight')
plt.xlabel('Log degree', fontsize=15)
plt.ylabel('Log frequency', fontsize=15)
plt.savefig('CESM_degree_distribution' + '.pdf', format='pdf')



plt.figure(figsize=(10.666666667, 6), dpi=200)
plt.xscale('log')
plt.yscale('log')
plt.scatter(np.linspace(1, len(hv), num=len(hv)), hv, c='#F67280', label='Hashimoto', marker='D')
plt.scatter(np.linspace(1, len(ev), num=len(ev)), ev, c='#355C7D', label='Eigen')
plt.xlim((0., len(ev)+10000))
plt.ylim((1e-24, 0.9))
plt.xlabel('Log rank', fontsize=15)
plt.ylabel('Log centrality', fontsize=15)
plt.legend()
plt.savefig('hashi_vs_eigen_cesm_abs_in' + '.pdf', format='pdf')