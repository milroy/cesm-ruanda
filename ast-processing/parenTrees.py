#!/bin/python

from metagraph_utils import *
from collections import defaultdict

class Node(object):
    """docstring for Node"""
    def __init__(self, depth, string, category):
        super(Node, self).__init__()
        self.category = category
        self.parent = None
        self.children = []
        self.depth = depth
        self.string = string
        self.inFromBelow = string
        self.isFunc = False
        self.isArray = False
        self.isParen = False
        self.defArgs = None
        self.isIface = False
        self.loc = None
 

class ParenTree(object):
    """docstring for ParenTree"""
    def __init__(self, string, topfunc, location, mgraph):
        super(ParenTree, self).__init__()
        self.location = copy.deepcopy(location)
        self.string = string
        self.nodes = defaultdict(list)

        for pcont in parenthetic_contents(self.string):
            node = Node(pcont[0], pcont[1])
            self.nodes[pcont[0]].append(node)

