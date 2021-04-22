import sys, os, re, copy
import fparser
from fparser import api
from fparser.api import Fortran2003
from fparser.readfortran import FortranStringReader
sys.path.insert(0, '/home/milroy/git/KGen/base/')
from kgen_extra import Intrinsic_Procedures
import networkx as nx
from subgraphs import metagraph
from metagraph_utils import *
import pdb

stmt = 'uflx(k)  = aleph(ie)%gimel(daled)+ zeta * foo((((aleph*wodka(bumble, psi)))) + cbmf * ( usrc + uplus - kp1 (  u0(kp1)  +   ssu0(kp1(10 + lambda(a), 1:100, epsilon=iota)) * ( ps0(k) - p0(kp1) + alpha(beta + max (gamma, nu) )) ) )) + aleph(nie)%gimel(daledd)'

#stmt = 'a = max(a, max(b, max( alpha, beta(min(gamma)))))'

#stmt = 'uflx(k)  = zeta * foo((((aleph*wodka(bumble, psi))))) + aleph(ie)%gimel(daled)'

reader = FortranStringReader(stmt)
assnStmt = Fortran2003.Assignment_Stmt(reader)

location = ['graph']
graph = metagraph('graph')

graph.functions['foo'] = {'args': ['one'], 'vars': ['a', 'b', 'c'], 'loc': ['graph', 'foomod'], 
                                          'output': 'foobar', 'isIface': False, 'isIntrinsic': False}

graph.functions['wodka'] = {'args': ['two', 'three'], 'vars': ['a', 'b', 'c'], 'loc': ['graph', 'wodkamod'], 
                                          'output': 'vodka', 'isIface': False, 'isIntrinsic': False}

graph.functions['u0'] = {'args': ['four'], 'vars': ['a', 'b', 'c'], 'loc': ['graph', 'u0mod'], 
                                          'output': 'u0out', 'isIface': False, 'isIntrinsic': False}

graph.functions['alpha'] = {'args': ['five'], 'vars': ['a', 'b', 'c'], 'loc': ['graph', 'alphamod'], 
                                          'output': 'alphaout', 'isIface': False, 'isIntrinsic': False}

graph.functions['max'] = {'args': None, 'vars': ['a', 'b', 'c'], 'loc': ['graph'], 
                                          'output': 'max', 'isIface': False, 'isIntrinsic': True}

lhs = assnStmt.items[0].tostr().split('(', 1)[0]
lhsSplit = lhs.split('(', 1)[0].strip().lower()
rhs = assnStmt.items[2].tostr().strip().lower()

lhsname = lhsSplit + '__' + location[-1]
lhsattrs = {'loc': location, 'cname': lhsSplit}
names, functions = stmt_to_fnnames(rhs)

names, functions = stmt_to_fnnames(rhs)
pVars = primitiveVars(graph, rhs, functions, location)

def getFuncAttrs(fname, location, mgraph):
    location = copy.deepcopy(location)
    isIface = False
    funcDef = None
    if fname in mgraph.functions:
        funcDef = mgraph.functions[fname]
    else:
        for idx, _ in enumerate(location):      
            keyLoc = tuple([fname]) + tuple(location[::-1][idx:][::-1])
            if keyLoc in mgraph.useMaps:
                mappedFname = mgraph.useMaps[keyLoc]
                funcDef = mgraph.functions[mappedFname]
    if funcDef:
        if not funcDef['isIface']:
            defArgs = funcDef['args']
            out = funcDef['output']
            loc = funcDef['loc']
            isIntrinsic = funcDef['isIntrinsic']
        else:
            defArgs = [funcDef['functions'][f]['args'] for f in funcDef['functions']]
            out = funcDef['functions'].itervalues().next()['output']
            loc = funcDef['functions'].itervalues().next()['loc']
            isIface = True
        return isIface, isIntrinsic, defArgs, out, loc
    else:
        raise KeyError(fname + ' not in CESM')

class Node(object):
    """docstring for Node"""
    def __init__(self, depth=None, string=None, order=None, location=location):
        super(Node, self).__init__()
        self.order = order
        self.parent = None
        self.parentFunction = None
        self.newNames = None
        self.newFunctions = None
        self.children = []
        self.depth = depth
        self.string = string
        self.strNoSpace = ''.join(string.split())
        self.inFromBelow = self.strNoSpace
        self.output = self.strNoSpace
        self.names, self.functions = stmt_to_fnnames(self.strNoSpace)
        self.isParen = False
        self.defArgs = None
        self.location = location
    def modifiedFunctions(self):
        self.newNames, self.newFunctions = stmt_to_fnnames(self.inFromBelow)
        self.output = self.inFromBelow

class ParenTree(object):
    """docstring for ParenTree"""
    def __init__(self, string, location, functionStruct):
        super(ParenTree, self).__init__()
        self.location = location#copy.deepcopy(location)
        self.string = string
        self.nodes = []
        self.levels = {}
        self.functions = functionStruct
        order = 1
        for lvl in parenthetic_contents(self.string):
            node = Node(depth=lvl[0], string=lvl[1], order=order, location=self.location)
            if lvl[0] not in self.levels:
                self.levels[lvl[0]] = [node]
            else:
                self.levels[lvl[0]].append(node)
            order += 1
        self.depth = len(self.levels)
        for i in reversed(xrange(1, self.depth)):
            for node in self.levels[i]:
                for nodeAbove in self.levels[i-1]:
                    if nodeAbove.order > node.order:
                        if '(' + node.strNoSpace + ')' in nodeAbove.strNoSpace:
                            node.parent = nodeAbove
                            nodeAbove.children.append(node)
                            for pf in nodeAbove.functions:
                                base = pf.split('(', 1)[0]
                                if pf == base + '(' + node.strNoSpace + ')':
                                    node.parentFunction = functionStruct[base]
                                    break
                            if not node.parentFunction:
                                # 
                                print 'paren'
                                node.parentFunction = ''
                                node.isParen = True
    def printNodes(self):
        for lvl, nodes in self.levels.iteritems():
            print "level: ", lvl
            for node in nodes:
                print "name: ", node.string
                print "children: ", [i.string for i in node.children]
                if node.parent:
                    print "parent: ", node.parent.string
                if node.parentFunction:
                    print "parent function: ", node.parentFunction.name
                elif node.parentFunction == '':
                    print "parent function: ''"




class Function(object):
    """docstring for Function"""
    def __init__(self, name, location, mgraph):
        super(Function, self).__init__()
        self.name = name
        isFunc = False
        isArray = False
        isParen = False
        defArgs = None
        isIface = False
        isIntrinsic = False
        if self.name in mgraph.functions:
            isFunc = True
            isIface, isIntrinsic, defArgs, out, loc = getFuncAttrs(self.name, location, mgraph)
        else:
            useMap = False
            for idx, _ in enumerate(location):
                keyLoc = tuple([self.name]) + tuple(location[::-1][idx:][::-1])
                if keyLoc in mgraph.useMaps:
                    target = mgraph.useMaps[keyLoc]
                    isIface, isIntrinsic, defArgs, out, loc = getFuncAttrs(target, location, mgraph)
                    isFunc = True
                    useMap = True
                    break
            if not useMap:
                loc = location
                if self.name == '' or self.name is not None:
                    out = self.name
                    isArray = True
                else:
                    out = None
                    isParen = True
        self.defArgs = defArgs
        self.output= out
        self.isFunc = isFunc
        self.isIface = isIface
        self.isArray = isArray
        self.isParen = isParen
        self.isIntrinsic = isIntrinsic
        self.location = loc


class FunctionTree(object):
    """docstring for FunctionTree"""
    def __init__(self, rootfunc, string, location, mgraph):
        super(FunctionTree, self).__init__()
        self.mgraph = mgraph
        self.string = string
        self.location = location#copy.deepcopy(location)
        names, self.functions = stmt_to_fnnames(self.string)
        self.funcStruct = {}
        for f in self.functions:
            fname = f.split('(', 1)[0]
            function = Function(fname, self.location, mgraph)
            self.funcStruct[fname] = function
        self.rootName = rootfunc
        self.root = Function(self.rootName, self.location, mgraph)
        self.funcStruct[self.rootName] = self.root
        self.argTree = ParenTree(string, self.location, self.funcStruct)
        # Set the root parenthetic's parent function to root function
        self.argTree.levels[0][0].parentFunction = self.root
        if len(self.argTree.levels[0]) > 1:
            raise RuntimeError("More than one root!")
    def argEdges(self, node, debug=False):
        argList = node.inFromBelow.split(',')
        newArgs = []
        for arg in argList:
            # Need this for precessed expressions,
            # e.g. foo(alpha*beta, gamma*epsilon)
            subnames, _ = stmt_to_fnnames(arg)
            newArgs.append(subnames)
        if not node.parentFunction.isArray:
            # Check for intrinsics
            if not node.parentFunction.isIntrinsic:
                if not node.parentFunction.isIface:
                    for defArg, inputArg in zip(node.parentFunction.defArgs, newArgs):
                        for subname in inputArg: 
                            if subname in node.newNames:
                                dargname = defArg + '__' + node.parentFunction.name
                                dargattrs = {'loc': node.parentFunction.location, 'cname': defArg}
                                snname = subname + '__' + node.location[-1]
                                snattrs = {'loc': node.location, 'cname': subname}
                                if debug:
                                    print "regular function: ", snname, dargname, snattrs, dargattrs
                                self.mgraph.addEdge(snname, dargname, attrs1=snattrs, attrs2=dargattrs)
                else:
                    # Map all argsets
                    for defArgSet in node.parentFunction.defArgs:
                        for defArg, inputArg in zip(defArgSet, newArgs):
                            for subname in inputArg: 
                                if subname in node.newNames:
                                    dargname = defArg + '__' + node.parentFunction.name
                                    dargattrs = {'loc': node.parentFunction.location, 'cname': defArg}
                                    snname = subname + '__' + node.location[-1]
                                    snattrs = {'loc': node.location, 'cname': subname}
                                    if debug:
                                        print "interface function: ", snname, dargname, snattrs, dargattrs
                                    self.mgraph.addEdge(snname, dargname, attrs1=snattrs, attrs2=dargattrs)
            else:
                # is intrinsic procedure
                for arg in newArgs:
                    for sn in arg: 
                        if sn in node.newNames:
                            #print darg, sn
                            funname = node.parentFunction.output + '__' + node.location[-1]
                            funattrs = {'loc': node.location, 'cname': node.parentFunction.output}
                            snname = sn + '__' + node.location[-1]
                            snattrs = {'loc': node.location, 'cname': sn}
                            if debug:
                                print "intrinsic procedure: ", snname, funname, snattrs, funattrs
                            self.mgraph.addEdge(snname, funname, attrs1=snattrs, attrs2=funattrs)
        else:
            # it's an array
            for arg in newArgs:
                for sn in arg: 
                    if sn in node.newNames:
                        #print shortestParent, sn
                        arrayname = node.parentFunction.output + '__' + node.location[-1]
                        arrayattrs = {'loc': node.location, 'cname': node.parentFunction.output}
                        snname = sn + '__' + node.location[-1]
                        snattrs = {'loc': node.location, 'cname': sn}
                        if debug:
                            print "array: ", snname, arrayname, snattrs, arrayattrs
                        self.mgraph.addEdge(snname, arrayname, attrs1=snattrs, attrs2=arrayattrs)
    def connectArgTree(self, debug=False):
        for lvl in reversed(xrange(1, self.argTree.depth)):
            for node in self.argTree.levels[lvl]:
                if node.parentFunction == '':
                    # eliminate one set of superfluous parentheses
                    node.output = node.inFromBelow
                    node.parent.inFromBelow = node.parent.inFromBelow.replace('(' + node.strNoSpace + ')', node.output, 1)
                elif node.parentFunction:
                    # node.connectFunctionArgs()
                    # node.parent.inFromBelow = node.output
                    node.modifiedFunctions()
                    argList = node.newNames
                    if node.parentFunction.defArgs and len(argList) != len(node.parentFunction.defArgs):
                        print argList, node.parentFunction.defArgs, node.parentFunction.name, node.string
                        #raise RuntimeError('input and definition arg lists different lengths')
                    self.argEdges(node, debug=debug)
                    # now rewrite string
                    node.parent.inFromBelow = node.parent.inFromBelow.replace(node.parentFunction.name + '(' + node.strNoSpace + ')', node.parentFunction.output, 1)
                    # and add edge
                    print node.parentFunction.name + '(' + node.strNoSpace + ')', node.parentFunction.output, node.parent.inFromBelow
        for node in self.argTree.levels[0]:
            node.modifiedFunctions()
            self.argEdges(node, debug=debug)
        return self.root.output


def connectfuncargs(parent, subargs, names, funcLike, location, mgraph):

    def processFArgs(defArgs, subargs, names, funcLike, floc, location, mgraph):
        floc = copy.deepcopy(floc)
        location = copy.deepcopy(location)
        defArgs = copy.deepcopy(defArgs)

        for darg, iarg in zip(defArgs, subargs):
            subnames, subfn = stmt_to_fnnames(iarg)

            # if subfunc, already added to fstdict, since doing DFS.
            if subfn:
                nodes = set([i.split('(', 1)[0] for i in subfn])

                for n in nodes:
                    suboutput = funcLike[n]['output']

                    if funcLike[n]['isFunc']:
                        _, _, _, subloc = getFuncAttrs(n, location, mgraph)

                    else:
                        subloc = location

                    #print darg, suboutput
                    dargname = darg + '__' + parent
                    dargattrs = {'loc': floc, 'cname': darg}
                    suboutname = suboutput + '__' + subloc[-1]
                    suboutattrs = {'loc': subloc, 'cname': suboutput}
                    #print suboutname, dargname, suboutattrs, dargattrs
                    mgraph.addEdge(suboutname, dargname, attrs1=suboutattrs, attrs2=dargattrs)

            else:
                for sn in subnames: 
                    if sn in names:
                        #print darg, sn
                        dargname = darg + '__' + parent
                        dargattrs = {'loc': floc, 'cname': darg}
                        snname = sn + '__' + location[-1]
                        snattrs = {'loc': location, 'cname': sn}
                        #print snname, dargname, snattrs, dargattrs
                        mgraph.addEdge(snname, dargname, attrs1=snattrs, attrs2=dargattrs)

        return None


    location = copy.deepcopy(location)
    isIface, defArgs, parout, floc = getFuncAttrs(parent, location, mgraph)

    #pdb.set_trace()

    # Need to check if is intrinsic function
    if defArgs is not None:

        if isIface:
            for dA in defArgs:
                processFArgs(dA, subargs, names, funcLike, floc, location, mgraph)

        else:
            processFArgs(defArgs, subargs, names, funcLike, floc, location, mgraph)

    else:

        for iarg in subargs:
            try:
                subnames, subfn = stmt_to_fnnames(iarg)
            except IndexError:
                #pdb.set_trace()
                pass


            # if subfunc, already added to fstdict, since doing DFS.
            if subfn:
                nodes = set([i.split('(', 1)[0] for i in subfn])

                for n in nodes:

                    try:
                        suboutput = funcLike[n]['output']

                        if funcLike[n]['isFunc']:
                             _, _, _, subloc = getFuncAttrs(n, location, mgraph)

                        else:
                            subloc = location

                    except KeyError:
                        suboutput = n
                        subloc = location
                        #pdb.set_trace()
                        pass


                    #print darg, suboutput
                    funname = parent + '__' + floc[-1]
                    funattrs = {'loc': floc, 'cname': parent}
                    suboutname = suboutput + '__' + subloc[-1]
                    suboutattrs = {'loc': subloc, 'cname': suboutput}
                    #print funname, suboutname
                    mgraph.addEdge(suboutname, funname, attrs1=suboutattrs, attrs2=funattrs)

            else:
                for sn in subnames: 
                    if sn in names:
                        #print darg, sn
                        funname = parent + '__' + floc[-1]
                        funattrs = {'loc': floc, 'cname': parent}
                        snname = sn + '__' + location[-1]
                        snattrs = {'loc': location, 'cname': sn}
                        #print funname, snname
                        mgraph.addEdge(snname, funname, attrs1=snattrs, attrs2=funattrs)

    return None





levels = {}
order = 1
for lvl in parenthetic_contents(self.string):
    if lvl[0] not in levels:
        levels[lvl[0]] = [(''.join(lvl[1].split()), lvl[1], order)]
    else:
        levels[lvl[0]].append([(''.join(lvl[1].split()), lvl[1], order)])
    order += 1

for i in reversed(xrange(1, len(dd))):
  for j in dd[i]:
    #print dd[i][j]
    for k in dd[i-1]:
      if k > j:
        if '(' + dd[i][j][0] + ')' in dd[i-1][k][0]:
          print i, j, k, dd[i][j][0], dd[i][j][1], dd[i-1][k][0], dd[i-1][k][1]



for i in reversed(xrange(1, len(dd))):
  for j in dd[i]:
    #print dd[i][j]
    for k in dd[i-1]:
      if k > j:
        for pa in parenthetic_contents(dd[i-1][k][0]):
          if dd[i][j][0] == pa[1]:
            print i, j, k, dd[i][j][0], pa, dd[i-1][k][0]




for lvl in reversed(xrange(pt.depth)):
  print lvl
  for node in pt.levels[lvl]:
    print 'string: ', node.string, node.order
    if node.parent:
      print 'parent: ', node.parent.string, node.parent.order
    print 'children'
    if node.children:
      for c in node.children:
        print c.string, c.order
