#!/bin/python

try:
    import graph_tool.all as g
    
except ImportError:
    import networkx as g

import sys
import re
import copy
import csv
#from suffix_trees import STree
from collections import deque
from collections import defaultdict
from collections import OrderedDict
from metagraph_utils import *
import pdb
from tqdm import tqdm
import traceback
import operator
from orderedset import OrderedSet
from itertools import tee
from itertools import combinations
from itertools import islice
import numpy as np
from math import sqrt
from BTrees.OOBTree import OOBTree

def pairwise(iterable):
    # From Python Docs
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


class metagraph(object):
    """This class creates and mutates NetworkX or
    graph-tool graphs.  It is designed to create 
    graphs of graphs to handle Fortran scoping in 
    the NCAR CESM."""

    def __init__(self, name=None):
        super(metagraph, self).__init__()

        if 'graph_tool' in sys.modules:
            """ TODO: replace with methods for graph_tool
            """ 
            self.graph = g.Graph(name=name)
            self.nodes = self.graph.nodes
            self.addnode = self.graph.add_node
            self.hasnode = self.graph.has_node
            self.addedge = self.graph.add_edge
            self.flatten = g.utils.flatten
            pass

        elif 'networkx' in sys.modules:
            self.graph = g.DiGraph(name=name)
            self.name = self.graph.name
            self.node = self.graph.node
            self.nodes = self.graph.nodes
            self.numNodes = self.graph.number_of_nodes
            self._addnode = self.graph.add_node
            self._addnodesfrom = self.graph.add_nodes_from
            self._remnodesfrom = self.graph.remove_nodes_from
            self._hasnode = self.graph.has_node
            self._hasedge = self.graph.has_edge
            self._addedge = self.graph.add_edge
            self.edges = self.graph.edges
            self._addedgesfrom = self.graph.add_edges_from
            self.numEdges = self.graph.number_of_edges
            self.flatten = g.utils.flatten
            self.subgraph = self.graph.subgraph
            self.outedges = self.graph.out_edges
            self.inedges = self.graph.in_edges

        self.cnames = bidict()
        self.modules = {}
        self.args = []
        self.calls = []
        self.functions = {}
        self.subroutines = {}
        self.funcLike = []
        self.uses = []
        self.useMaps = {}
        self.interfaces = []
        self.outputs = set()
        self.inputs = set()

        removeChars = ['+', '/', '-', '*', '==', '<', '>']
        self.removeWords = ['.and.', '.or.']
        self.rx = '[' + re.escape(''.join(removeChars)) + ']'

        self.camVars = set(['trefht', 'icimr', 'burdenseasalt',
        'bc_a1_srf', 'ts', 'numliq', 'psl', 'h2o2_srf', 'burdenbc',
        'cldtot', 'freqs', 'freqr', 'tgcldiwp', 'cldice', 'freql',
        'so4_a2_srf', 't', 'freqi', 'flut', 'cldmed', 'cdnumc',
        'qrefht', 'numice', 'tmq', 'so4_a3_srf', 'num_a1_srf',
        'cloud', 'tgcldlwp', 'burdendust', 'icwmr', 'fsnt', 'aodvis',
        'soa_a1_srf', 'taux', 'ncl_a3_srf', 'z3', 'aoddust1',
        'aoddust3', 'fsdsc', 'dst_a1_srf', 'cldhgh', 'wgustd', 'fsds',
        'flnsc', 'qrl', 'precsc', 'pom_a1_srf', 'qrs', 'num_a2_srf',
        'flntc', 'ps', 'lhflx', 'arel', 'flns', 'flnt', 'fsntoac',
        'arei', 'anrain', 'u10', 'fsnsc', 'so4_a1_srf', 'ncl_a1_srf',
        'taugwx', 'taugwy', 'soa_a2_srf', 'flds', 'aqsnow', 'so2_srf',
        'relhum', 'aqrain', 'awnc', 'uu', 'fsns', 'v', 'ansnow',
        'tauy', 'dtv', 'fsntc', 'pblh', 'num_a3_srf', 'fsntoa',
        'cldlow', 'precsl', 'dtcond', 'vq', 'wsub', 'vt', 'vu', 'vv',
        'precl', 'precc', 'burdensoa', 'awni', 'fice', 'omegat',
        'iwc', 'ncl_a2_srf', 'flutc', 'burdenso4', 'dst_a3_srf',
        'omega', 'snowhlnd', 'shflx', 'qflx', 'ccn3', 'burdenpom',
        'q', 'h2so4_srf', 'u', 'cldliq'])

        self.outMaps = {'dtv': 'dtk__vertical_diffusion_tend',
        'lhflx': 'cam_in%lhf__diag_surf', 'qflx':
        'cam_in%cflx__diag_surf', 'qrefht': 'cam_in%qref__diag_surf',
        'shflx': 'cam_in%shf__diag_surf', 'snowhlnd':
        'cam_in%snowhland__diag_surf', 'taugwx': 'tau0x__gw_tend',
        'taugwy': 'tau0y__gw_tend', 'taux': 'cam_in%wsx__diag_surf',
        'tauy': 'cam_in%wsy__diag_surf', 'trefht':
        'cam_in%tref__diag_surf', 'u10': 'cam_in%u10__diag_surf',
        'wgustd': 'wpert__vertical_diffusion_tend', 'ts':
        'cam_in%ts__diag_surf', 'flds':
        'cam_out%flwds__radiation_tend'}

        with open('module-file-lists/atmfiles.csv', 'r') as f:
            reader = csv.reader(f, delimiter=',')
            self.atmFiles = set(list(reader)[0])

        with open('module-file-lists/csm_sharefiles.csv', 'r') as f:
            reader = csv.reader(f, delimiter=',')
            self.csm_shareFiles = set(list(reader)[0])

        with open('module-file-lists/dead_sharefiles.csv', 'r') as f:
            reader = csv.reader(f, delimiter=',')
            self.dead_shareFiles = set(list(reader)[0])

        with open('module-file-lists/drvfiles.csv', 'r') as f:
            reader = csv.reader(f, delimiter=',')
            self.drvFiles = set(list(reader)[0])

        with open('module-file-lists/glcfiles.csv', 'r') as f:
            reader = csv.reader(f, delimiter=',')
            self.glcFiles = set(list(reader)[0])

        with open('module-file-lists/icefiles.csv', 'r') as f:
            reader = csv.reader(f, delimiter=',')
            self.iceFiles = set(list(reader)[0])

        with open('module-file-lists/lndfiles.csv', 'r') as f:
            reader = csv.reader(f, delimiter=',')
            self.lndFiles = set(list(reader)[0])

        with open('module-file-lists/ocnfiles.csv', 'r') as f:
            reader = csv.reader(f, delimiter=',')
            self.ocnFiles = set(list(reader)[0])

        with open('module-file-lists/roffiles.csv', 'r') as f:
            reader = csv.reader(f, delimiter=',')
            self.rofFiles = set(list(reader)[0])

        with open('module-file-lists/utilsfiles.csv', 'r') as f:
            reader = csv.reader(f, delimiter=',')
            self.utilsFiles = set(list(reader)[0])

        with open('module-file-lists/wavfiles.csv', 'r') as f:
            reader = csv.reader(f, delimiter=',')
            self.wavFiles = set(list(reader)[0])

        with open('module-file-lists/onlycam.csv', 'r') as f:
            reader = csv.reader(f, delimiter=',')
            self.camFiles = set(list(reader)[0])

    """The following function is taken from 
    https://stackoverflow.com/questions/10823877/what-is-the-fastest-way-to-flatten-arbitrarily-nested-lists-in-python
    """
    def __flatten(self, container):
        for i in container:
            if isinstance(i, (list, tuple)):
                for j in self.__flatten(i):
                    yield j
            else:
                yield i

    def _parse_eq(self, string, removeX, removeW):
        tmp = re.sub(removeX, ' ', string.lower()).strip().split()
        resultWords  = [word for word in tmp if word not in removeW]
        return resultWords

    def _parenContents(self, string):
        """
        From stackoverflow question 4284991:
        Generate parenthesized contents in string as pairs (level, contents).
        """
        stack = []
        for i, c in enumerate(string):
            if c == '(':
                stack.append(i)
            elif c == ')' and stack:
                start = stack.pop()
                yield len(stack), string[start + 1: i]

    def subgraphs(self):
        sgList = self.subGraphs.keys()
        for s in self.subGraphs:
            sgSGL = len(self.subGraphs[s].subgraphs())
            if sgSGL > 0:
                sgList.append(self.subGraphs[s].subgraphs())

        return sgList   

    def addNode(self, node, attrs=None):
        attributes = copy.deepcopy(attrs)
        baseline = attributes['lines']['line']
        if 'line' in attributes['lines']:
            # this happens for usemaps defined in the module "header"
            if not isinstance(baseline, OOBTree):
                attributes['lines']['line'] = OOBTree({baseline: 0})

        if attributes['cname']:
            if '%' in node:
                base = node.split('__')[0].split('%')[-1]
                attributes['cname'] = base
                self.cnames[node] = base

            else:
                self.cnames[node] = attributes['cname']

        self._addnode(node, attrs=attributes)

        return None

    def addEdge(self, n1, n2, attrs1=None, attrs2=None, edgeLine=None):
        attributes1 = copy.deepcopy(attrs1)
        attributes2 = copy.deepcopy(attrs2)

        if '%' in n1:
            if '(' in n1:
                n1 = re.sub(r'\([^)]*\)', '', n1)
                if ')' in n1:
                    n1 = n1.replace(')', '')

        if '%' in n2:
            if '(' in n2:
                n2 = re.sub(r'\([^)]*\)', '', n2)
                if ')' in n2:
                    n2 = n2.replace(')', '')

        if not self._hasnode(n1):
            self.addNode(n1, attrs=attributes1)

        else:
            if 'lines' in attributes1:
                if isinstance(attributes1['lines']['line'], OOBTree):
                    self.node[n1]['attrs']['lines']['line'].update(attributes1['lines']['line'])
                elif isinstance(attributes1['lines']['line'], int):
                    self.node[n1]['attrs']['lines']['line'].update({attributes1['lines']['line']: 0})

        if not self._hasnode(n2):
            self.addNode(n2, attrs=attributes2)

        else:
            if 'lines' in attributes2:
                if isinstance(attributes2['lines']['line'], OOBTree):
                    self.node[n2]['attrs']['lines']['line'].update(attributes2['lines']['line'])
                elif isinstance(attributes2['lines']['line'], int):
                    self.node[n2]['attrs']['lines']['line'].update({attributes2['lines']['line']: 0})

        if not self._hasedge(n1, n2):
            if edgeLine:
                if isinstance(edgeLine, OOBTree):
                    self._addedge(n1, n2, lines=edgeLine)
                elif isinstance(edgeLine, int):
                    self._addedge(n1, n2, lines=OOBTree({edgeLine: 0}))
                else:
                    pdb.set_trace()

            else:
                pdb.set_trace()

        else:
            if edgeLine:
                if isinstance(edgeLine, OOBTree):
                    self.edges[n1, n2]['lines'].update(edgeLine)
                elif isinstance(edgeLine, int):
                    self.edges[n1, n2]['lines'].update({edgeLine: 0})
                else:
                    pdb.set_trace()

            else:
                pdb.set_trace()

        return None


    def addNodesFrom(self, ebunch, attrbunch=None):
        cebunch = copy.deepcopy(ebunch)
        cattrbunch = copy.deepcopy(attrbunch)
        flatbunch = self.flatten(cebunch)
        if not len(flatbunch) == len(cattrbunch):
            print "Warning: assuming nodes already have attributes"

            self._addnodesfrom(flatbunch)

        for node, attrs in zip(flatbunch, cattrbunch):
            self.addNode(node, attrs=attrs)

        return None

    def addEdgesFrom(self, ebunch, attrbunch):
        cebunch = copy.deepcopy(ebunch)
        cattrbunch = copy.deepcopy(attrbunch)        
        self.addNodesFrom(cebunch, cattrbunch)
        self._addedgesfrom(cebunch)

        return None

    def getParent(self, node):
        if node in self.subGraphs:
            return self.name
        else:
            parent = None
            for sg in self.subGraphs:
                if self.subGraphs[sg].hasNode(node):
                    parent = self.subGraphs[sg].getParent(node)

            return parent


    def hasLike(self, subString):
        where = []
        nodes = self.nodes()
        nodes = filter(None, nodes)

        return [var for var in nodes if subString in var]


    def hasPath(self, source=None, target=None):
        return g.has_path(self.graph, source=source, target=target)

    def shortestPath(self, source=None, target=None, weight=None):
        return g.shortest_path(self.graph, source=source, target=target, weight=weight)

    def processUsemaps(self):
        for use in self.uses:

            useLocation = use[2]
            for item in use[1]['renamed']:

                # Fortran use rename is confusing: targetName => sourceName
                if item[1] in self.functions or item[1] in self.subroutines:
                    locKey = tuple([item[0]]) + tuple(useLocation)
                    self.useMaps[locKey] = item[1]

        return None

    def processUsedVars(self):
        # need to do this last to make sure all vars exist in graph
        for use in self.uses:
            useLocation = use[2][-1]
            originModule = use[0]

            if use[1]['renamed']:
                # Fortran use rename is confusing: targetName => sourceName
                for item in use[1]['renamed']:
                    sourceName = item[1]
                    targetName = item[0]

                    # check if subprogram
                    if sourceName not in self.functions and sourceName not in self.subroutines:

                        if originModule in self.modules:

                            if sourceName in self.modules[originModule]['vars'] | self.modules[originModule]['public'] | self.modules[originModule]['provides']:
                                smlocation = self.modules[originModule]['location']
                                matchedSources = [(sourceName + '__' + originModule, {'location': smlocation, 'cname': sourceName, 'lines': {'line': OOBTree({1: 0})}})]

                                for subprog in self.modules[originModule]['subprogs']:
                                    suffix = '__' + subprog

                                    if sourceName + suffix in self.nodes: 

                                        newSourceName = sourceName + suffix

                                        if subprog in self.subroutines:
                                            if 'location' in self.subroutines[subprog]:
                                                location = self.subroutines[subprog]['location']

                                            elif self.subroutines[subprog]['isIface']:
                                                subiface = self.subroutines[subprog]['subroutines'].itervalues().next()
                                                location = subiface['location']

                                            else:
                                                location = None

                                        elif subprog in self.functions:
                                            if 'location' in self.functions[subprog]:
                                                location = self.functions[subprog]['location']

                                            elif self.functions[subprog]['isIface']:
                                                subiface = self.functions[subprog]['functions'].itervalues().next()
                                                location = subiface['location']

                                            else:
                                                location = None

                                        else:
                                            location = None

                                        if location:
                                            lines = self.node[newSourceName]['attrs']['lines']
                                            sourceattrs = {'location': location, 'cname': sourceName, 'lines': lines}
                                            matchedSources.append((newSourceName, sourceattrs))

                                matchedTargets = []
                                if useLocation in self.modules:
                                    # Need to look for all variables matching string in module subprograms
                                    for subprog in self.modules[useLocation]['subprogs']:
                                        suffixT = '__' + subprog

                                        if targetName + suffixT in self.nodes:

                                            newTargetName = targetName + suffixT

                                            if subprog in self.subroutines:
                                                if 'location' in self.subroutines[subprog]:
                                                    location = self.subroutines[subprog]['location']

                                                elif self.subroutines[subprog]['isIface']:
                                                    subiface = self.subroutines[subprog]['subroutines'].itervalues().next()
                                                    location = subiface['location']

                                                else:
                                                    location = None

                                            elif subprog in self.functions:
                                                if 'location' in self.functions[subprog]:
                                                    location = self.functions[subprog]['location']

                                                elif self.functions[subprog]['isIface']:
                                                    subiface = self.functions[subprog]['functions'].itervalues().next()
                                                    location = subiface['location']

                                                else:
                                                    location = None

                                            else:
                                                location = None

                                            if location:
                                                lines = self.node[newTargetName]['attrs']['lines']
                                                targetattrs = {'location': location, 'cname': targetName, 'lines': lines}
                                                matchedTargets.append((newTargetName, targetattrs))

                                elif useLocation in self.subroutines:
                                    suffixT = '__' + useLocation

                                    if targetName + suffixT in self.nodes:

                                        newTargetName = targetName + suffixT
                                        lines = self.node[newTargetName]['attrs']['lines']
                                        targetattrs = {'location': self.subroutines[useLocation]['location'], 'cname': targetName, 'lines': lines}
                                        matchedTargets.append((newTargetName, targetattrs))

                                elif useLocation in self.functions:
                                    suffixT = '__' + useLocation

                                    if targetName + suffixT in self.nodes:

                                        newTargetName = targetName + suffixT
                                        lines = self.node[newTargetName]['attrs']['lines']
                                        targetattrs = {'location': self.functions[useLocation]['location'], 'cname': targetName, 'lines': lines}
                                        matchedTargets.append((newTargetName, targetattrs))

                                if matchedTargets and matchedSources:

                                    for mt in matchedTargets:
                                        for ms in matchedSources:

                                            self.addEdge(ms[0], mt[0], attrs1=ms[1], attrs2=mt[1], edgeLine=mt[1]['lines']['line'])

                            else:
                                print 'not public ', item, use

                        else:
                            print originModule, ' not in self modules'

            if use[1]['same']:
                for name in use[1]['same']:

                    sourceName = name
                    targetName = name

                    # check if subprogram
                    if name not in self.functions and name not in self.subroutines:

                        if originModule in self.modules:

                            if sourceName in self.modules[originModule]['vars'] | self.modules[originModule]['public'] | self.modules[originModule]['provides']:
                                smlocation = self.modules[originModule]['location']
                                matchedSources = [(sourceName + '__' + originModule, {'location': smlocation, 'cname': sourceName, 'lines': {'line': OOBTree({1: 0})}})]
                                for subprog in self.modules[originModule]['subprogs']:
                                    suffix = '__' + subprog

                                    if sourceName + suffix in self.nodes: 

                                        newSourceName = sourceName + suffix

                                        if subprog in self.subroutines:
                                            if 'location' in self.subroutines[subprog]:
                                                location = self.subroutines[subprog]['location']

                                            elif self.subroutines[subprog]['isIface']:
                                                subiface = self.subroutines[subprog]['subroutines'].itervalues().next()
                                                location = subiface['location']

                                            else:
                                                location = None

                                        elif subprog in self.functions:
                                            if 'location' in self.functions[subprog]:
                                                location = self.functions[subprog]['location']

                                            elif self.functions[subprog]['isIface']:
                                                subiface = self.functions[subprog]['functions'].itervalues().next()
                                                location = subiface['location']

                                            else:
                                                location = None

                                        else:
                                            location = None

                                        if location:
                                            lines = self.node[newSourceName]['attrs']['lines']
                                            sourceattrs = {'location': location, 'cname': sourceName, 'lines': lines}
                                            matchedSources.append((newSourceName, sourceattrs))

                                matchedTargets = []
                                if useLocation in self.modules:
                                    # Need to look for all variables matching string in module subprograms
                                    for subprog in self.modules[useLocation]['subprogs']:
                                        suffixT = '__' + subprog

                                        if targetName + suffixT in self.nodes:

                                            newTargetName = targetName + suffixT

                                            if subprog in self.subroutines:
                                                if 'location' in self.subroutines[subprog]:
                                                    location = self.subroutines[subprog]['location']

                                                elif self.subroutines[subprog]['isIface']:
                                                    subiface = self.subroutines[subprog]['subroutines'].itervalues().next()
                                                    location = subiface['location']

                                                else:
                                                    location = None

                                            elif subprog in self.functions:
                                                if 'location' in self.functions[subprog]:
                                                    location = self.functions[subprog]['location']

                                                elif self.functions[subprog]['isIface']:
                                                    subiface = self.functions[subprog]['functions'].itervalues().next()
                                                    location = subiface['location']

                                                else:
                                                    location = None

                                            else:
                                                location = None

                                            if location:
                                                lines = self.node[newTargetName]['attrs']['lines']
                                                targetattrs = {'location': location, 'cname': targetName, 'lines': lines}
                                                matchedTargets.append((newTargetName, targetattrs))

                                elif useLocation in self.subroutines:
                                    suffixT = '__' + useLocation

                                    if targetName + suffixT in self.nodes:

                                        newTargetName = targetName + suffixT
                                        lines = self.node[newTargetName]['attrs']['lines']
                                        targetattrs = {'location': self.subroutines[useLocation]['location'], 'cname': targetName, 'lines': lines}
                                        matchedTargets.append((newTargetName, targetattrs))

                                elif useLocation in self.functions:
                                    suffixT = '__' + useLocation

                                    if targetName + suffixT in self.nodes:

                                        newTargetName = targetName + suffixT
                                        lines = self.node[newTargetName]['attrs']['lines']
                                        targetattrs = {'location': self.functions[useLocation]['location'], 'cname': targetName, 'lines': lines}
                                        matchedTargets.append((newTargetName, targetattrs))

                                if matchedTargets and matchedSources:

                                    for mt in matchedTargets:
                                        for ms in matchedSources:

                                            self.addEdge(ms[0], mt[0], attrs1=ms[1], attrs2=mt[1], edgeLine=mt[1]['lines']['line'])

                            else:
                                print 'not public ', item, use

                        else:
                            print originModule, ' not in self modules'

            elif use[1]['all']:

                if originModule in self.modules:

                    for var in self.modules[originModule]['vars'] | self.modules[originModule]['public'] | self.modules[originModule]['provides']:
                        sourceName = var
                        targetName = var

                        # check if subprogram
                        if var not in self.functions and var not in self.subroutines:

                            if originModule in self.modules:

                                if sourceName in self.modules[originModule]['vars'] | self.modules[originModule]['public'] | self.modules[originModule]['provides']:
                                    smlocation = self.modules[originModule]['location']
                                    matchedSources = [(sourceName + '__' + originModule, {'location': smlocation, 'cname': sourceName, 'lines': {'line': OOBTree({1: 0})}})]
                                    for subprog in self.modules[originModule]['subprogs']:
                                        suffix = '__' + subprog

                                        if sourceName + suffix in self.nodes: 

                                            newSourceName = sourceName + suffix

                                            if subprog in self.subroutines:
                                                if 'location' in self.subroutines[subprog]:
                                                    location = self.subroutines[subprog]['location']

                                                elif self.subroutines[subprog]['isIface']:
                                                    subiface = self.subroutines[subprog]['subroutines'].itervalues().next()
                                                    location = subiface['location']

                                                else:
                                                    location = None

                                            elif subprog in self.functions:
                                                if 'location' in self.functions[subprog]:
                                                    location = self.functions[subprog]['location']

                                                elif self.functions[subprog]['isIface']:
                                                    subiface = self.functions[subprog]['functions'].itervalues().next()
                                                    location = subiface['location']

                                                else:
                                                    location = None

                                            else:
                                                location = None

                                            if location:
                                                lines = self.node[newSourceName]['attrs']['lines']
                                                sourceattrs = {'location': location, 'cname': sourceName, 'lines': lines}
                                                matchedSources.append((newSourceName, sourceattrs))

                                    matchedTargets = []
                                    if useLocation in self.modules:
                                        # Need to look for all variables matching string in module subprograms
                                        for subprog in self.modules[useLocation]['subprogs']:
                                            suffixT = '__' + subprog

                                            if targetName + suffixT in self.nodes:

                                                newTargetName = targetName + suffixT

                                                if subprog in self.subroutines:
                                                    if 'location' in self.subroutines[subprog]:
                                                        location = self.subroutines[subprog]['location']

                                                    elif self.subroutines[subprog]['isIface']:
                                                        subiface = self.subroutines[subprog]['subroutines'].itervalues().next()
                                                        location = subiface['location']

                                                    else:
                                                        location = None

                                                elif subprog in self.functions:
                                                    if 'location' in self.functions[subprog]:
                                                        location = self.functions[subprog]['location']

                                                    elif self.functions[subprog]['isIface']:
                                                        subiface = self.functions[subprog]['functions'].itervalues().next()
                                                        location = subiface['location']

                                                    else:
                                                        location = None

                                                else:
                                                    location = None

                                                if location:
                                                    lines = self.node[newTargetName]['attrs']['lines']
                                                    targetattrs = {'location': location, 'cname': targetName, 'lines': lines}
                                                    matchedTargets.append((newTargetName, targetattrs))

                                    elif useLocation in self.subroutines:
                                        suffixT = '__' + useLocation

                                        if targetName + suffixT in self.nodes:

                                            newTargetName = targetName + suffixT
                                            lines = self.node[newTargetName]['attrs']['lines']
                                            targetattrs = {'location': self.subroutines[useLocation]['location'], 'cname': targetName, 'lines': lines}
                                            matchedTargets.append((newTargetName, targetattrs))

                                    elif useLocation in self.functions:
                                        suffixT = '__' + useLocation

                                        if targetName + suffixT in self.nodes:

                                            newTargetName = targetName + suffixT
                                            lines = self.node[newTargetName]['attrs']['lines']
                                            targetattrs = {'location': self.functions[useLocation]['location'], 'cname': targetName, 'lines': lines}
                                            matchedTargets.append((newTargetName, targetattrs))

                                    if matchedTargets and matchedSources:

                                        for mt in matchedTargets:
                                            for ms in matchedSources:

                                                self.addEdge(ms[0], mt[0], attrs1=ms[1], attrs2=mt[1], edgeLine=mt[1]['lines']['line'])

                                else:
                                    print 'not public ', item, use

                            else:
                                print originModule, ' not in self modules'

        return None


    def processIfaces(self):
        for iface in self.interfaces:

            procedure = {}
            procedure['isIface'] = True
            for subprog in iface[1]:

                if subprog in self.subroutines:
                    if not 'subroutines' in procedure:
                        procedure['subroutines'] = {subprog: self.subroutines[subprog]}

                    else:
                        procedure['subroutines'][subprog] = self.subroutines[subprog]

                elif subprog in self.functions:
                    if not 'functions' in procedure:
                        procedure['functions'] = {subprog: self.functions[subprog]}

                    else:
                        procedure['functions'][subprog] = self.functions[subprog]

                else:
                    print ValueError(subprog + ' not in CESM')
                    pass

            # Must be subroutine or function.
            if subprog in self.subroutines:
                try:
                    self.subroutines[iface[0]] = procedure

                except:
                    raise KeyError(iface[0] + ' already in CESM subroutines')

            if subprog in self.functions:
                try:
                    self.functions[iface[0]] = procedure

                except:
                    raise KeyError(iface[0] + ' already in CESM functions')

        return None

    def connectCallArgs(self):

        errs = [0, []]
        for statement in tqdm(self.calls, unit='subrcalls'):

            if not statement[3]['processed']:

                subname = statement[0]
                callArgs = statement[1]
                location = statement[2]

                subroutine = None

                try:

                    if subname in self.subroutines:
                        subroutine = self.subroutines[subname]

                    else:

                        for idx, _ in enumerate(location['subprogs']):
                        
                            keyLoc = tuple([subname]) + tuple(location['subprogs'][::-1][idx:][::-1])
                            if keyLoc in self.useMaps:
                                mappedSname = self.useMaps[keyLoc]
                                subroutine = self.subroutines[mappedSname]

                    if subroutine:

                        if subroutine['isIface']:

                            try:

                                for k, subsubr in subroutine['subroutines'].iteritems():
                                    subLoc = subsubr['location']
                                    for carg, defarg in zip(callArgs, subsubr['args']):
                                        processSubArgs(self, carg, defarg, location, subLoc)

                            except KeyError:
                                errs[0] += 1
                                errs[1].append(statement)
                                pass

                        else:
                            subLoc = subroutine['location']
                            for carg, defarg in zip(callArgs, subroutine['args']):
                                processSubArgs(self, carg, defarg, location, subLoc)


                    else:
                        errs[0] += 1
                        errs[1].append(statement)
                        #print traceback.format_exc()
                        #print exc
                        #pdb.set_trace()

                except Exception as e:
                    errs[0] += 1
                    errs[1].append(statement)

        return errs


    def connectFunctions(self):

        errs = [0, []]
        for statement in tqdm(self.funcLike, unit='fcalls'):

            try:

                location = statement[1]
                stmt = statement[0]

                if isinstance(stmt, str):
                    LHS = stmt.split('=', 1)[0].strip().lower()
                    if '%' not in LHS:
                        lhsSplit = stmt.split('=', 1)[0].strip().lower().split('(', 1)[0]

                    else:
                        lhsSplit = re.sub(r'\([^)]*\)', '', LHS)
                        if ')' in lhsSplit:
                            lhsSplit = lhsSplit.replace(')', '')

                    processedLHS = ''.join(lhsSplit.split())

                    rhs = stmt.split('=', 1)[1].strip().lower()

                else:
                    LHS = stmt.items[0].tostr()
                    if '%' not in LHS:
                        lhsSplit = LHS.split('(', 1)[0].strip().lower()

                    else:
                        lhsSplit = re.sub(r'\([^)]*\)', '', LHS.strip().lower())
                        if ')' in lhsSplit:
                            lhsSplit = lhsSplit.replace(')', '')

                    processedLHS = ''.join(lhsSplit.split())

                    rhs = stmt.items[2].tostr().strip().lower()

                lhsname = processedLHS + '__' + location['subprogs'][-1]
                lhsattrs = {'location': location['subprogs'], 'cname': processedLHS, 'lines': location['lines']}
                names, functions = stmt_to_fnnames(rhs)

                rootVs = rootVars(self, rhs, functions, location)

                if ' ' in lhsname:
                    print stmt, lhsname

                for rv in rootVs:
                    self.addEdge(rv[0], lhsname, attrs1=rv[1], attrs2=lhsattrs, edgeLine=lhsattrs['lines']['line'])

                if functions:
                    deepStmts(self, rhs, location, lhsname, lhsattrs)

            except Exception as e:
                errs[0] += 1
                errs[1].append(statement)
                #pdb.set_trace()
                pass

        return errs


    def isectPaths(self, name, lastn=0, outputs=False, cnames=True):
        pathSet = set()
        nodeSet = set()

        if cnames:
            for cn in self.cnames.inverse[name]:
                if not outputs:
                    for s, v in self.shortestPath(target=cn).iteritems():
                        for source, target in pairwise(v):
                            pathSet.add((source, target))
                            nodeSet.add(source)

                else:
                    if cn in self.outputs:
                        for s, v in self.shortestPath(target=cn).iteritems():
                            for source, target in pairwise(v):
                                pathSet.add((source, target))
                                nodeSet.add(source)
                nodeSet.add(cn)

        else:
            for n in self.hasLike(name):
                if not outputs:
                    for s, v in self.shortestPath(target=n).iteritems():
                        for source, target in pairwise(v):
                            pathSet.add((source, target))
                            nodeSet.add(source)

                else:
                    if n in self.outputs:
                        for s, v in self.shortestPath(target=n).iteritems():
                            for source, target in pairwise(v):
                                pathSet.add((source, target))
                                nodeSet.add(source)
                nodeSet.add(n)

        return pathSet, nodeSet


    def isectLocs(self, name):
        hlSet = set()

        for nlike in self.hasLike(name):
            for ie in self.inedges(nlike):
                hlSet.update(self.node[ie[0]]['attrs']['location'])
                    
        return hlSet

    def getSubgraph(self, name):
        nodes = [n for n in self.nodes() if name in self.node[n]['attrs']['location'][0]]
        
        return g.subgraph(self, nodes)

    def __resolveNames(self, name, inedges=True, weight=None):

        hlSet = set()
        for nlike in self.hasLike(name):
            for p in self.cnames.inverse[name]:

                if self.hasPath(nlike, p) or self.hasPath(p, nlike):
                   
                    if inedges:
                        for ie in self.inedges(nlike):
                            hlSet.add(ie[0])
                            if weight:
                                self.edges[ie[0], nlike]['weight'] = weight

                    else:
                        for oe in self.edges(nlike):
                            hlSet.add(oe[1])
                            if weight:
                                self.edges[nlike, oe[1]]['weight'] = weight

                    hlSet.add(nlike)

        return hlSet 

    def hashimoto(self, sgraph, inout=False, outCentrality=False, weight=None, max_iter=50, tol=0):
        import scipy as sp

        if len(sgraph) == 0:
            raise g.NetworkXPointlessConcept('cannot compute centrality for the'
                                              ' null graph')

        def hashimoto_matrix(G, create_using=None):
            """
            Modified to return the Hashimoto non-backtracking matrix.

            Return the line graph L of the (multi)digraph G.

            Edges in G appear as nodes in L, represented as tuples of the form (u,v)
            or (u,v,key) if G is a multidigraph. A node in L corresponding to the edge
            (u,v) is connected to every node corresponding to an edge (v,w).

            Parameters
            ----------
            G : digraph
                A directed graph or directed multigraph.
            create_using : None
                A digraph instance used to populate the line graph.

            """
            if create_using is None:
                L = G.fresh_copy()
            else:
                L = create_using

            # Create a graph specific edge function.
            get_edges = g.generators.line._edge_func(G)

            for from_node in get_edges():
                # from_node is: (u,v) or (u,v,key)
                L.add_node(from_node)
                for to_node in get_edges(from_node[1]):
                    # Modification here: make sure there aren't reciprocal edges.
                    if to_node != (from_node[1], from_node[0]):
                        L.add_edge(from_node, to_node)
            return L

        if not outCentrality:

            hashimoto_mat = hashimoto_matrix(sgraph)

            L = g.to_scipy_sparse_matrix(hashimoto_mat, nodelist=list(hashimoto_mat), dtype=float)
            eigenvalue, eigenvector = sp.sparse.linalg.eigs(L.T, k=1, which='LR', maxiter=50, tol=0)

            # Need absolute value, since the norm can be complex. Must investigate.
            largest = np.abs(eigenvector.flatten().real)
            # norm = largest.sum()
            # norml = largest/norm
            # now find the nodes' centralities
            ncenters = defaultdict(list)
            for idx, n in enumerate(hashimoto_mat):
                # line graph nodes are edge tuples of original graph
                # want inedges by definition of centrality
                ncenters[n[1]].append(idx)


        if inout or outCentrality:
            rsg = sgraph.reverse()

            rhashimoto_mat = hashimoto_matrix(rsg)

            rL = g.to_scipy_sparse_matrix(rhashimoto_mat, nodelist=list(rhashimoto_mat), dtype=float)
            reigenvalue, reigenvector = sp.sparse.linalg.eigs(rL.T, k=1, which='LR', maxiter=50, tol=0)

            # Need absolute value, since the norm can be complex. Must investigate.
            rlargest = np.abs(reigenvector.flatten().real)
            # rnorm = rlargest.sum()
            # rnorml = rlargest/rnorm
            # now find the nodes' centralities
            rncenters = defaultdict(list)
            for idx, n in enumerate(rhashimoto_mat):
                # line graph nodes are edge tuples of original graph
                # want inedges by definition of centrality
                rncenters[n[1]].append(idx)


        centrality = {}
        if not outCentrality and not inout:
            tmp = []
            for n in sorted(ncenters):
                tmp.append(np.sum(largest[ncenters[n]]))

            normsign = np.sign(largest.sum())
            tmpnorm = normsign*np.array(tmp, dtype=np.float64)/np.linalg.norm(tmp)
            for idx, n in enumerate(sorted(ncenters)):
                centrality[n] = tmpnorm[idx]

        elif outCentrality and not inout:
            tmp = []
            for n in sorted(rncenters):
                tmp.append(np.sum(rlargest[rncenters[n]]))

            rnormsign = np.sign(rlargest.sum())
            tmpnorm = rnormsign*np.array(tmp, dtype=np.float64)/np.linalg.norm(tmp)
            for idx, n in enumerate(sorted(rncenters)):
                centrality[n] = tmpnorm[idx]

        else:
            tmpin = tmpout = []
            # need to combine unique keys
            combined = set(ncenters.keys() + rncenters.keys())

            for n in sorted(combined):
                tmpin.append(np.sum(largest[ncenters[n]]))
                tmpout.append(np.sum(rlargest[rncenters[n]]))

            normsign = np.sign(largest.sum())
            rnormsign = np.sign(rlargest.sum())
            tmpinNorm = normsign*np.array(tmpin, dtype=np.float64)/np.linalg.norm(tmpin)
            tmpoutNorm = rnormsign*np.array(tmpout, dtype=np.float64)/np.linalg.norm(tmpout)
            # since each centrality is normed to 1, norm sum to 1 by Cauchy-Schwartz:
            sumNorm = (tmpinNorm + tmpoutNorm)/sqrt(2. + np.inner(tmpinNorm, tmpoutNorm) + np.inner(tmpoutNorm, tmpinNorm))
            for idx, n in enumerate(sorted(combined)):
                centrality[n] = sumNorm[idx]

            
        return centrality

    def girvan_newman(self, inducedSgraph, nComms=1, hashi=False, outCentrality=False, inout=False):
        
        gn = g.algorithms.community.centrality.girvan_newman(inducedSgraph)
        communities = None
        for comms in islice(gn, nComms):
            communities = tuple(sorted(c) for c in comms)

        if communities:
            comGraphs = []
            for c in communities:

                if len(c) > 1:

                    try:
                        csg = g.subgraph(self.graph, c)

                        if hashi:
                            eigens = self.hashimoto(csg, inout=inout, outCentrality=outCentrality)

                        else:
                            if outCentrality:
                                eigens = g.eigenvector_centrality_numpy(csg.reverse())

                            elif inout:
                                ieigenNodes = g.eigenvector_centrality_numpy(csg)
                                oeigenNodes = g.eigenvector_centrality_numpy(csg.reverse())
                                tmpin = tmpout = []
                                for n in csg:
                                    tmpin.append(ieigenNodes[n])
                                    tmpout.append(oeigenNodes[n])

                                tmpinNorm = np.array(tmpin, dtype=np.float64)/np.linalg.norm(tmpin)
                                tmpoutNorm = np.array(tmpout, dtype=np.float64)/np.linalg.norm(tmpout)
                                # since each centrality is normed to 1, norm to 1 by Cauchy-Schwartz:
                                sumNorm = (tmpinNorm + tmpoutNorm)/sqrt(2. + np.inner(tmpinNorm, tmpoutNorm) + np.inner(tmpoutNorm, tmpinNorm))

                                eigens = {}
                                for idx, n in enumerate(csg):
                                    eigens[n] = sumNorm[idx]

                            else:
                                eigens = g.eigenvector_centrality_numpy(csg)

                        sorted_eigens = sorted(eigens.items(), key=operator.itemgetter(1), reverse=True)
                        tmp = {'subgraph': csg, 'centralities': sorted_eigens}
                        comGraphs.append(tmp)

                    except:
                        pass

            return comGraphs

    def eigenCentrality(self, nameList, num=-1, nhubs=10, lastn=0, topLocs=10000000, inedges=True, \
                        hubs=False, paths=True, notVars=False, allCESM=False, outputs=False, \
                        cnames=True, overlap=False, models=None, quotient=False, qindex=0, \
                        both=False, qweights=False, hashi=False, inout=False, girvan=False, \
                        nComms=1, outCentrality=False):

        if models:
            mdict = {'atm': self.atmFiles, 'csm_share': self.csm_shareFiles, 'dead_share': self.dead_shareFiles, \
                     'drv': self.drvFiles, 'glc': self.glcFiles, 'ice': self.iceFiles, 'lnd': self.lndFiles, \
                     'ocn': self.ocnFiles, 'rof': self.rofFiles, 'utils': self.utilsFiles, \
                     'wav': self.wavFiles, 'cam': self.camFiles}

            mlist = []
            for m in models:
                mlist.append(mdict[m])

            modelMods = set.union(*mlist)

        if isinstance(nameList, dict):
            weights = True

        else:
            nameList = set(nameList)
            weights = False

        if not paths:
            if weights:
                inedgeLocs = [self.__resolveNames(k, weight=v) for k, v in nameList.iteritems()]

            else:
                inedgeLocs = [self.__resolveNames(k) for k in nameList]

            ninedgeLocs = set()
            if notVars:
                for k in self.camVars:
                    if k not in nameList and k not in set(['v', 'ps', 't', 'ts', 'q', 'u']):
                        ninedgeLocs.update(self.__resolveNames(k))

        else:
            if weights:
                inedgeLocs = [self.isectPaths(k, outputs=outputs, cnames=cnames)[0] for k, v in nameList.iteritems()]

                if notVars:
                    ninedgeLocs = {}
                    for k in self.camVars:
                        if k not in nameList and k not in set(['v', 'ps', 't', 'ts', 'q', 'u']):
                            ninedgeLocs[k] = self.isectPaths(k, outputs=outputs, cnames=cnames)[0]

            else:
                inedgeLocs = [self.isectPaths(i, outputs=outputs, cnames=cnames)[0] for i in nameList]

                if notVars:
                    ninedgeLocs = set()
                    for k in self.camVars:
                        if k not in nameList and k not in set(['v', 'ps', 't', 'ts', 'q', 'u']):
                            ninedgeLocs.update(self.isectPaths(k, outputs=outputs, cnames=cnames)[0])



        if not notVars and not allCESM:
            if paths:
                allNodes = set()
                for a, b in set.union(*inedgeLocs):
                    allNodes.add(a)
                    allNodes.add(b)

            else:
                allNodes = set.union(*inedgeLocs)

        elif notVars and not allCESM:
            if not paths:
                posNodes = set.union(*inedgeLocs)
                allNodes = posNodes - ninedgeLocs

            elif paths and weights:
                allNodes = set()

                if overlap:
                    for a, b in combinations(inedgeLocs, 2):
                        for ovrlp in a & b:
                            self.edges[ovrlp[0], ovrlp[1]]['weight'] *= 2.

                for a, b in set.union(*inedgeLocs):
                    allNodes.add(a)
                    allNodes.add(b)

                for diffpath in inedgeLocs:
                    for var, notPath in ninedgeLocs.iteritems():
                        for isect in diffpath & notPath:
                            self.edges[isect[0], isect[1]]['weight'] /= 2.

            else:
                posNodes = set.union(*inedgeLocs)
                allEdges = posNodes - ninedgeLocs

                allNodes = set()
                for a, b in allEdges:
                    allNodes.add(a)
                    allNodes.add(b)

        elif allCESM:
            allNodes = set(self.nodes)


        if models:
            if paths:
                for a, b in set.union(*inedgeLocs):
                    src = self.node[a]['attrs']['location'][0]
                    dest = self.node[b]['attrs']['location'][0]

                    if src not in modelMods and dest not in modelMods:
                        allNodes.discard(a)
                        allNodes.discard(b)

            else:
                for n in allNodes:
                    if self.node[n]['attrs']['location'][0] not in modelMods:
                        allNodes.remove(n)

        if not quotient:
            inducedSgraph = g.subgraph(self.graph, allNodes)


        else:

            eqRel = lambda u, v: (self.node[u]['attrs']['location'][qindex] == self.node[v]['attrs']['location'][qindex])

            sg = g.subgraph(self.graph, allNodes)

            if not qweights:
                inducedSgraph = g.quotient_graph(sg, eqRel)

            else:
                def edgeRel(U, V):
                    weight = 0
                    for node in U:
                        for vert in V:
                            if self.graph.has_edge(node, vert) or self.graph.has_edge(vert, node):
                                weight += 1

                    return {'weight': weight}


                inducedSgraph = g.quotient_graph(sg, eqRel, edge_data=edgeRel)


        print "subgraph nodes:", inducedSgraph.number_of_nodes(), " edges:", inducedSgraph.number_of_edges()

        if inedges:
            if weights:
                if not hashi:
                    eigenNodes = g.eigenvector_centrality_numpy(inducedSgraph, weight='weight')

                else:
                    eigenNodes = self.hashimoto(inducedSgraph, inout=inout, outCentrality=outCentrality, weight='weight')

            else:
                if not hashi:
                    if inout:
                        ieigenNodes = g.eigenvector_centrality_numpy(inducedSgraph)
                        oeigenNodes = g.eigenvector_centrality_numpy(inducedSgraph.reverse())
                        eigenNodes = {}
                        for n in inducedSgraph:
                            eigenNodes[n] = abs(ieigenNodes[n]) + abs(oeigenNodes[n])

                    elif outCentrality:
                        eigenNodes = g.eigenvector_centrality_numpy(inducedSgraph.reverse())

                    else:
                        eigenNodes = g.eigenvector_centrality_numpy(inducedSgraph)

                else:
                    eigenNodes = self.hashimoto(inducedSgraph, inout=inout, outCentrality=outCentrality)

        else:
            if weights:
                eigenNodes = g.eigenvector_centrality_numpy(inducedSgraph.reverse(), weight='weight')

            else:
                eigenNodes = g.eigenvector_centrality_numpy(inducedSgraph.reverse())

        sorted_eigen = sorted(eigenNodes.items(), key=operator.itemgetter(1), reverse=True)
        # save the first centrality for function return
        firstEigen = copy.deepcopy(sorted_eigen)

        importantLocs = OrderedDict()
        importantNodes = set()
        top2 = OrderedSet()

        if not quotient:

            for i in sorted_eigen[:num]:
                if self.name not in i[0]:
                    top2.update(self.node[i[0]]['attrs']['location'])
                    if i[0] not in importantNodes:
                        importantNodes.add(i[0])
                        key = tuple(self.node[i[0]]['attrs']['location'])
                        if key not in importantLocs:
                            importantLocs[key] = i[1]                    

            if hubs:
                allNodes.remove(sorted_eigen[0][0])
                inducedSgraph = g.subgraph(self.graph, allNodes)
                print "subgraph nodes:", inducedSgraph.number_of_nodes(), " edges:", inducedSgraph.number_of_edges()

                for _ in xrange(1, nhubs):

                    try:

                        if inedges:
                            if weights:
                                eigenNodes = g.eigenvector_centrality_numpy(inducedSgraph, weight='weight')

                            else:
                                eigenNodes = g.eigenvector_centrality_numpy(inducedSgraph)

                        else:
                            if weights:
                                eigenNodes = g.eigenvector_centrality_numpy(inducedSgraph.reverse(), weight='weight')

                            else:
                                eigenNodes = g.eigenvector_centrality_numpy(inducedSgraph.reverse())

                    except:
                        break

                    sorted_eigen = sorted(eigenNodes.items(), key=operator.itemgetter(1), reverse=True)
                    for i in sorted_eigen[:num]:
                        if self.name not in i[0]:
                            top2.update(self.node[i[0]]['attrs']['location'])
                            if i[0] not in importantNodes:
                                importantNodes.add(i[0])
                                key = tuple(self.node[i[0]]['attrs']['location'])
                                if key not in importantLocs:
                                    importantLocs[key] = i[1]

                    allNodes.remove(sorted_eigen[0][0])
                    inducedSgraph = g.subgraph(self.graph, allNodes)
                    print "subgraph nodes:", inducedSgraph.number_of_nodes(), " edges:", inducedSgraph.number_of_edges()

            # Reset edge weights
            if paths and weights and notVars:
                print "Resetting weights"
                for path in inedgeLocs:
                    for pathEdge in path:
                        self.edges[pathEdge[0], pathEdge[1]]['weight'] = 1


        else: # Need to operate on quotient graphs

            for i in sorted_eigen[:num]:
                if isinstance(i[0], frozenset):
                    nodeName = iter(i[0]).next()

                else:
                    nodeName = i[0]

                try:

                    if self.name not in nodeName:
                        if 'origin' in self.node[nodeName]['attrs']:
                            top2.update(self.node[nodeName]['attrs']['origin'])
                            if nodeName not in importantNodes:
                                importantNodes.add(nodeName)
                                key = tuple(self.node[nodeName]['attrs']['origin'])
                                if key not in importantLocs:
                                    importantLocs[key] = i[1]
                        else:
                            top2.update(self.node[nodeName]['attrs']['location'])
                            if nodeName not in importantNodes:
                                importantNodes.add(nodeName)
                                key = tuple(self.node[nodeName]['attrs']['location'])
                                if key not in importantLocs:
                                    importantLocs[key] = i[1]

                except:
                    pdb.set_trace()
              

            if hubs:

                inducedSgraph.remove_node(sorted_eigen[0][0])
                print "subgraph nodes:", inducedSgraph.number_of_nodes(), " edges:", inducedSgraph.number_of_edges()

                for _ in xrange(1, nhubs):

                    try:

                        if inedges:
                            if weights:
                                eigenNodes = g.eigenvector_centrality_numpy(inducedSgraph, weight='weight')

                            else:
                                eigenNodes = g.eigenvector_centrality_numpy(inducedSgraph)

                        else:
                            if weights:
                                eigenNodes = g.eigenvector_centrality_numpy(inducedSgraph.reverse(), weight='weight')

                            else:
                                eigenNodes = g.eigenvector_centrality_numpy(inducedSgraph.reverse())

                    except:
                        break

                    sorted_eigen = sorted(eigenNodes.items(), key=operator.itemgetter(1), reverse=True)
                    for i in sorted_eigen[:num]:
                        if isinstance(i[0], frozenset):
                            nodeName = iter(i[0]).next()

                        else:
                            nodeName = i[0]

                        if self.name not in nodeName:
                            if 'origin' in self.node[nodeName]['attrs']:
                                top2.update(self.node[nodeName]['attrs']['origin'])
                                if nodeName not in importantNodes:
                                    importantNodes.add(nodeName)
                                    key = tuple(self.node[nodeName]['attrs']['origin'])
                                    if key not in importantLocs:
                                        importantLocs[key] = i[1]
                            else:
                                top2.update(self.node[nodeName]['attrs']['location'])
                                if nodeName not in importantNodes:
                                    importantNodes.add(nodeName)
                                    key = tuple(self.node[nodeName]['attrs']['location'])
                                    if key not in importantLocs:
                                        importantLocs[key] = i[1]

                    inducedSgraph.remove_node(sorted_eigen[0][0])
                    print "subgraph nodes:", inducedSgraph.number_of_nodes(), " edges:", inducedSgraph.number_of_edges()


        top = OrderedSet()
        for k in sorted(importantLocs.items(), key=operator.itemgetter(1), reverse=True):
            if len(top) < topLocs:
                top.add(k[0][0])
            else:
                break

        unsrt = set()
        for m in importantLocs:
            unsrt.add(m[0])

        if girvan:
            comGraphs = self.girvan_newman(inducedSgraph, hashi=hashi, nComms=nComms, outCentrality=outCentrality, inout=inout)

            if comGraphs:
                return top, top2, unsrt, importantLocs, allNodes, firstEigen, comGraphs

            else:
                return top, top2, unsrt, importantLocs, allNodes, firstEigen

        else:
            return top, top2, unsrt, importantLocs, allNodes, firstEigen


    def backwardSlice(self, node):
        backslices = set()
        for path in self.shortestPath(target=node).items():
            spstack = {}
            currentSP = None
            #print path[1]
            rpath = list(reversed(path[1]))
            target = rpath[0]
            targetLine = self.node[target]['attrs']['lines']['line'].maxKey(sys.maxsize)
            backslice = [(target, targetLine)]
            idx = 2
            lenrpath = len(rpath)
            for source in rpath[1:]:
                #print source, target
                ssubProg = source.split('__')[-1]
                if ssubProg == currentSP:
                    assignLines = set.intersection(*[set(self.edges[source, target]['lines']), \
                        set(self.node[source]['attrs']['lines']['line']), \
                        set(self.node[target]['attrs']['lines']['line'])])
                    assignOOB = OOBTree({i: 0 for i in assignLines})
                    try:
                        nextLine = assignOOB.maxKey(spstack[ssubProg])
                        spstack[ssubProg] = nextLine
                        backslice.append((source, nextLine))
                        target = source
                        currentSP = ssubProg
                    except ValueError:
                        backslices.add(tuple(reversed(backslice)))
                        break
                else:
                    currentSP = ssubProg
                    if ssubProg not in spstack:
                        nextLine = self.node[source]['attrs']['lines']['line'].maxKey(sys.maxsize)
                        spstack[ssubProg] = nextLine
                        backslice.append((source, nextLine))
                        target = source
                    else:
                        try:
                            nextLine = self.node[source]['attrs']['lines']['line'].maxKey(spstack[ssubProg] -1)
                            spstack[ssubProg] = nextLine
                            backslice.append((source, nextLine))
                            target = source
                        except ValueError:
                            backslices.add(tuple(reversed(backslice)))
                            break

                if idx == lenrpath:
                    backslices.add(tuple(reversed(backslice)))

                else:
                    idx += 1

        return backslices


    def forwardSlice(self, node):
        backslices = set()
        for path in self.shortestPath(source=node).items():
            spstack = {}
            currentSP = None
            #print path[1]
            rpath = list((path[1]))
            source = rpath[0]
            forwardslice = [source]
            idx = 2
            lenrpath = len(rpath)
            for target in rpath[1:]:
                #print source, target
                ssubProg = target.split('__')[-1]
                if ssubProg == currentSP:
                    assignLines = set.intersection(*[set(self.edges[source, target]['lines']), \
                        set(self.node[source]['attrs']['lines']['line']), \
                        set(self.node[target]['attrs']['lines']['line'])])
                    assignOOB = OOBTree({i: 0 for i in assignLines})
                    try:
                        spstack[ssubProg] = assignOOB.minKey(spstack[ssubProg])
                        forwardslice.append(target)
                        source = target
                        currentSP = ssubProg
                    except ValueError:
                        backslices.add(tuple(forwardslice))
                        break
                else:
                    currentSP = ssubProg
                    if ssubProg not in spstack:
                        spstack[ssubProg] = self.node[target]['attrs']['lines']['line'].minKey(0)
                        forwardslice.append(target)
                        source = target
                    else:
                        try:
                            spstack[ssubProg] = self.node[source]['attrs']['lines']['line'].minKey(spstack[ssubProg] +1)
                            forwardslice.append(target)
                            source = target
                        except ValueError:
                            backslices.add(tuple(forwardslice))
                            break

                if idx == lenrpath:
                    backslices.add(tuple(forwardslice))

                else:
                    idx += 1

        return backslices