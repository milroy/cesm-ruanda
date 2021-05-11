#!/bin/python

import collections
import copy
import unicodedata
import fparser
from fparser import api
from fparser import Fortran2003
from fparser.readfortran import FortranStringReader
import pdb
import sys
import re
import psutil, os

# Try this for deep composite calls
sys.setrecursionlimit(2000)

class bidict(dict):
    """
    https://stackoverflow.com/questions/3318625/efficient-bidirectional-hash-table-in-python
    """
    def __init__(self, *args, **kwargs):
        super(bidict, self).__init__(*args, **kwargs)
        self.inverse = {}
        for key, value in self.iteritems():
            self.inverse.setdefault(value,[]).append(key)

    def __setitem__(self, key, value):
        if key in self:
            self.inverse[self[key]].remove(key) 
        super(bidict, self).__setitem__(key, value)
        self.inverse.setdefault(value,[]).append(key)

    def __delitem__(self, key):
        self.inverse.setdefault(self[key],[]).remove(key)
        if self[key] in self.inverse and not self.inverse[self[key]]: 
            del self.inverse[self[key]]
        super(bidict, self).__delitem__(key)


class Node(object):
    """
    Node in a tree of parenthetic statements
    """
    def __init__(self, depth=None, string=None, order=None, location=None):
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


class Function(object):
    """
    Node in a tree of Functions and argument tree
    """
    def __init__(self, name, location, mgraph):
        super(Function, self).__init__()
        self.name = name
        isFunc = False
        isArray = False
        isParen = False
        defArgs = None
        isIface = False
        isIntrinsic = False
        fTarget = None
        useMap = False

        if self.name in mgraph.functions:
            isFunc = True
            isIface, isIntrinsic, defArgs, out, loc = getFuncAttrs(self.name, location['subprogs'], mgraph)

        else:
            for idx, _ in enumerate(location['subprogs']):
                keyLoc = tuple([self.name]) + tuple(location['subprogs'][::-1][idx:][::-1])
                if keyLoc in mgraph.useMaps:
                    target = mgraph.useMaps[keyLoc]
                    isIface, isIntrinsic, defArgs, out, loc = getFuncAttrs(target, location['subprogs'], mgraph)
                    isFunc = True
                    useMap = True
                    fTarget = target
                    self.name = target
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
        self.isUsed = useMap
        self.fTarget = fTarget


class ParenTree(object):
    """
    Class for methods and attributes of deep parenthetic statements
    """
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
                                node.parentFunction = ''
                                node.isParen = True

    def printNodes(self):
        # Print for debugging
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


class FunctionTree(object):
    """
    Wrapper class for Functions, Nodes, and ParenTrees
    """
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
        if len(self.argTree.levels[0]) > 1:
            # ok if derived type
            if '%' in self.string:
                for rootNode in self.argTree.levels[0]:
                    rootNode.parentFunction = self.root
            else:
                raise RuntimeError("More than one root!")
        else:
            self.argTree.levels[0][0].parentFunction = self.root
        

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
                                lines = self.mgraph.node[dargname]['attrs']['lines']
                                dargattrs = {'location': node.parentFunction.location, 'cname': defArg, 'lines': lines}
                                snname = subname + '__' + node.location['subprogs'][-1]
                                snattrs = {'location': node.location['subprogs'], 'cname': subname, 'lines': node.location['lines']}

                                if debug:
                                    print "regular function: ", snname, dargname, snattrs, dargattrs

                                self.mgraph.addEdge(snname, dargname, attrs1=snattrs, attrs2=dargattrs, edgeLine=snattrs['lines']['line'])
                else:
                    # Map all argsets
                    for defArgSet in node.parentFunction.defArgs:
                        for defArg, inputArg in zip(defArgSet, newArgs):
                            for subname in inputArg: 
                                if subname in node.newNames:
                                    # Need to find which interface subprogram has dargname
                                    if node.parentFunction.name in self.mgraph.functions:
                                        for child in self.mgraph.functions[node.parentFunction.name]['functions']:
                                            if defArg in self.mgraph.functions[node.parentFunction.name]['functions'][child]['vars']:
                                                dargname = defArg + '__' + child
                                                lines = self.mgraph.node[dargname]['attrs']['lines']
                                                dargattrs = {'location': self.mgraph.functions[node.parentFunction.name]['functions'][child]['location'], 'cname': defArg, 'lines': lines}
                                                break
                                    elif node.parentFunction.name in self.mgraph.subroutines:
                                        b = False
                                        for child in self.mgraph.subroutines[node.parentFunction.name]['subroutines']:
                                            for arg in self.mgraph.subroutines[node.parentFunction.name]['subroutines'][child]['args']:
                                                if defArg == arg[0]:
                                                    dargname = defArg + '__' + child
                                                    lines = self.mgraph.node[dargname]['attrs']['lines']
                                                    dargattrs = {'location': self.mgraph.subroutines[node.parentFunction.name]['subroutines'][child]['location'], 'cname': defArg, 'lines': lines}
                                                    b = True
                                                    break
                                            if b:
                                                break
                                    else:
                                        # Need to add code to find the used interface function (this will take a while).
                                        dargname = defArg + '__' + node.parentFunction.name
                                        dargattrs = {'location': node.parentFunction.location, 'cname': defArg, 'lines': node.location['lines']}

                                    snname = subname + '__' + node.location['subprogs'][-1]
                                    snattrs = {'location': node.location['subprogs'], 'cname': subname, 'lines': node.location['lines']}

                                    if debug:
                                        print "interface function: ", snname, dargname, snattrs, dargattrs

                                    self.mgraph.addEdge(snname, dargname, attrs1=snattrs, attrs2=dargattrs, edgeLine=snattrs['lines']['line'])
            else:
                # is intrinsic procedure
                for arg in newArgs:
                    for sn in arg: 
                        if sn in node.newNames:
                            #print darg, sn
                            funname = node.parentFunction.output + '__' + str(node.location['lines']['line']) + '__' + node.location['subprogs'][-1]
                            funattrs = {'location': node.location['subprogs'], 'cname': node.parentFunction.output, 'lines': node.location['lines']}
                            snname = sn + '__' + node.location['subprogs'][-1]
                            snattrs = {'location': node.location['subprogs'], 'cname': sn, 'lines': node.location['lines']}

                            if debug:
                                print "intrinsic procedure: ", snname, funname, snattrs, funattrs

                            self.mgraph.addEdge(snname, funname, attrs1=snattrs, attrs2=funattrs, edgeLine=snattrs['lines']['line'])
        else:
            # it's an array
            pass


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
                    # if node.parentFunction.defArgs and len(argList) != len(node.parentFunction.defArgs):
                    #     print argList, node.parentFunction.defArgs, node.parentFunction.name, node.string
                    #     #raise RuntimeError('input and definition arg lists different lengths')
                    self.argEdges(node, debug=debug)
                    # now rewrite string
                    node.parent.inFromBelow = node.parent.inFromBelow.replace(node.parentFunction.name + '(' + node.strNoSpace + ')', node.parentFunction.output, 1)
                    # and add edge
                    # print node.parentFunction.name + '(' + node.strNoSpace + ')', node.parentFunction.output, node.parent.inFromBelow

        for node in self.argTree.levels[0]:
            node.modifiedFunctions()
            self.argEdges(node, debug=debug)


class OrderedSet(collections.MutableSet):
    """
    https://code.activestate.com/recipes/576694/
    """

    def __init__(self, iterable=None):
        self.end = end = [] 
        end += [None, end, end]         # sentinel node for doubly linked list
        self.map = {}                   # key --> [key, prev, next]
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def add(self, key):
        if key not in self.map:
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]

    def discard(self, key):
        if key in self.map:        
            key, prev, next = self.map.pop(key)
            prev[2] = next
            next[1] = prev

    def __iter__(self):
        end = self.end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = self.end[1][0] if last else self.end[2][0]
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)


def get_names(node, bag, f2003_type, depth):
    if isinstance(node, f2003_type):
        bag.append(''.join(node.string.split()))


def traverse(node, func, bag, f2003_type, subnode='items', prerun=True, depth=0):
    ret = None
    if prerun and func is not None:
        ret = func(node, bag, f2003_type, depth)
        if ret is not None: return ret

    if node and hasattr(node, subnode) and getattr(node, subnode) is not None:
        for child in getattr(node, subnode):
            ret = traverse(child, func, bag, f2003_type, subnode=subnode, prerun=prerun, depth=depth+1)

    if not prerun and func is not None:
        ret = func(node, bag, f2003_type, depth)
        if ret is not None: return ret

    return ret


def parenthetic_contents(string):
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

            yield (len(stack), string[start + 1: i])


def stmt_to_fnnames(instmt):
    names = []
    functions = []
    arrays = []
    zerofunc = []
    derivedTypes = []

    if isinstance(instmt, str):

        if not instmt:
            return set([]), set([])

        else:

            # dummy = 'INVALID = INVALID(' + instmt.replace('!', '') + ')'
            dummy = 'INVALID = INVALID(' + instmt + ')'
            reader = FortranStringReader(dummy)
            stmt = Fortran2003.Assignment_Stmt(reader)
            traverse(stmt, get_names, derivedTypes, fparser.Fortran2003.Data_Ref)
            if derivedTypes:
                dtnew = []
                for dt in derivedTypes:
                    if '(' in dt:
                        # remove indices, since we don't map them anyway
                        dtnospace = re.sub(r'\([^)]*\)', '', ''.join(dt.split()))
                        if ')' in dtnospace:
                            dtnospace = dtnospace.replace(')', '')

                    else:
                        dtnospace = ''.join(dt.split())
                    
                    dtnew.append(dtnospace)

                string = ''.join(stmt.tostr().split())
                dtmap = {}
                num = 0
                for dt, dtn in zip(derivedTypes, dtnew):
                    dtnospace = ''.join(dt.split())
                    # We must order the replacements so we don't replace with the wrong string
                    string = string.replace(dtnospace, 'REMOVETHISSTRING' + '_' + str(num) + '_', 1)
                    dtmap[num] = dtn
                    num += 1

                derivedTypes = dtnew

                reader = FortranStringReader(string)
                stmt = Fortran2003.Assignment_Stmt(reader)

            traverse(stmt, get_names, names, fparser.Fortran2003.Name)
            traverse(stmt, get_names, functions, fparser.Fortran2003.Part_Ref)
            traverse(stmt, get_names, zerofunc, fparser.Fortran2003.Structure_Constructor)
            traverse(stmt, get_names, arrays, fparser.Fortran2003.Array_Section)
            functions.extend(zerofunc)
            functions.extend(arrays)
            functions = [f for f in functions if 'INVALID' not in f]
            if derivedTypes:
                newNames = []
                for n in names:
                    if 'REMOVETHISSTRING' in n:
                        idx = int(n.split('_')[1])
                        newNames.append(dtmap[idx])

                    else:
                        newNames.append(n)

                names = newNames

                newFunctions = []
                for f in functions:
                    newf = f                       
                    if 'REMOVETHISSTRING' in f:
                        matches = re.findall('REMOVETHISSTRING_[0-9]+_', f)
                        for m in matches:
                            idx = int(m.split('_')[1])
                            newf = newf.replace(m, dtmap[idx], 1)

                    newFunctions.append(newf)

                functions = newFunctions

            names = OrderedSet(names[2:])
            functions = set(functions)
    else:

        traverse(instmt, get_names, derivedTypes, fparser.Fortran2003.Data_Ref)
        if derivedTypes:
            dtnew = []
            for dt in derivedTypes:
                if '(' in dt:
                    # remove indices, since we don't map them anyway
                    dtnospace = re.sub(r'\([^)]*\)', '', ''.join(dt.split()))
                    if ')' in dtnospace:
                        dtnospace = dtnospace.replace(')', '')

                else:
                    dtnospace = ''.join(dt.split())
                
                dtnew.append(dtnospace)

            string = ''.join(instmt.tostr().split())
            dtmap = {}
            num = 0
            for dt, dtn in zip(derivedTypes, dtnew):
                dtnospace = ''.join(dt.split())
                # We must order the replacements so we don't replace with the wrong string
                string = string.replace(dtnospace, 'REMOVETHISSTRING' + '_' + str(num) + '_', 1)
                dtmap[num] = dtn
                num += 1

            derivedTypes = dtnew

            reader = FortranStringReader(string)
            instmt = Fortran2003.Assignment_Stmt(reader)

        traverse(instmt, get_names, names, fparser.Fortran2003.Name)
        traverse(instmt, get_names, functions, fparser.Fortran2003.Part_Ref)
        traverse(instmt, get_names, zerofunc, fparser.Fortran2003.Structure_Constructor)
        traverse(instmt, get_names, arrays, fparser.Fortran2003.Array_Section)
        functions.extend(zerofunc)
        functions.extend(arrays)
        if derivedTypes:
            newNames = []
            for n in names:
                if 'REMOVETHISSTRING' in n:
                    idx = int(n.split('_')[1])
                    newNames.append(dtmap[idx])

                else:
                    newNames.append(n)

            names = newNames

            newFunctions = []
            for f in functions:
                newf = f
                if 'REMOVETHISSTRING' in f:
                    matches = re.findall('REMOVETHISSTRING_[0-9]+_', f)
                    for m in matches:
                        idx = int(m.split('_')[1])
                        newf = newf.replace(m, dtmap[idx], 1)

                newFunctions.append(newf)

            functions = newFunctions

        names = OrderedSet(names)
        functions = set(functions)

    return names, functions


def getFuncAttrs(fname, location, mgraph):
    #location = copy.deepcopy(location)
    isIface = False
    isIntrinsic = False

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
            loc = funcDef['location']
            isIntrinsic = funcDef['isIntrinsic']

        else:
            defArgs = [funcDef['functions'][f]['args'] for f in funcDef['functions']]
            out = funcDef['functions'].itervalues().next()['output']
            loc = funcDef['functions'].itervalues().next()['location']
            isIface = True

        return isIface, isIntrinsic, defArgs, out, loc

    else:
        raise KeyError(fname + ' not in CESM')


def connectOutput(functionTree, location, mgraph, targetname=None, targetattrs=None, isOut=True):

    # location = copy.deepcopy(location)
    # targetname = copy.deepcopy(targetname)
    # targetattrs = copy.deepcopy(targetattrs)
    # arguments = output.split(',')
    output = functionTree.root.output

    if len(output.split(',')) > 1:
        raise RuntimeError('output of functionStruct is a list')

    if functionTree.root.isFunc:
        if functionTree.root.isIntrinsic:
            rname = functionTree.rootName
            rootLoc = location['subprogs']
            outname = output + '__' + str(location['lines']['line']) + '__' + location['subprogs'][-1]
            outattrs = {'location': rootLoc, 'cname': outname, 'lines': location['lines']}

        else:
            rname = functionTree.rootName
            rootLoc = functionTree.root.location

            outname = output + '__' + rootLoc[-1]
            outattrs = {'location': rootLoc, 'cname': outname}
        
        if isOut:
            mgraph.addEdge(outname, targetname, attrs1=outattrs, attrs2=targetattrs, edgeLine=targetattrs['lines']['line'])
            #mgraph.addEdge(targetname, outname, attrs1=lhsattrs, attrs2=outattrs)

        elif 'direction' in targetattrs:
            if targetattrs['direction'] == 'in':
                mgraph.addEdge(outname, targetname, attrs1=outattrs, attrs2=targetattrs, edgeLine=outattrs['lines']['line'])

            elif targetattrs['direction'] == 'out':
                mgraph.addEdge(targetname, outname, attrs1=targetattrs, attrs2=outattrs, edgeLine=outattrs['lines']['line'])

            # is 'inout'
            else:
                mgraph.addEdge(outname, targetname, attrs1=outattrs, attrs2=targetattrs, edgeLine=outattrs['lines']['line'])
                mgraph.addEdge(targetname, outname, attrs1=targetattrs, attrs2=outattrs, edgeLine=outattrs['lines']['line'])

        else:
            mgraph.addEdge(targetname, outname, attrs1=targetattrs, attrs2=outattrs, edgeLine=outattrs['lines']['line'])
            #mgraph.addEdge(outname, targetname, attrs1=outattrs, attrs2=targetattrs)

        return None

    elif functionTree.root.isArray:

        outname = output + '__' + location['subprogs'][-1]
        outattrs = {'location': location['subprogs'], 'cname': output, 'lines': location['lines']}

        if isOut:
            mgraph.addEdge(outname, targetname, attrs1=outattrs, attrs2=targetattrs, edgeLine=outattrs['lines']['line'])
            #mgraph.addEdge(targetname, outname, attrs1=targetattrs, attrs2=outattrs)

        elif 'direction' in targetattrs:
            #print targetattrs['direction']
            if targetattrs['direction'] == 'in':
                mgraph.addEdge(outname, targetname, attrs1=outattrs, attrs2=targetattrs, edgeLine=outattrs['lines']['line'])

            elif targetattrs['direction'] == 'out':
                mgraph.addEdge(targetname, outname, attrs1=targetattrs, attrs2=outattrs, edgeLine=outattrs['lines']['line'])

            # is 'inout'
            else:
                mgraph.addEdge(outname, targetname, attrs1=outattrs, attrs2=targetattrs, edgeLine=outattrs['lines']['line'])
                mgraph.addEdge(targetname, outname, attrs1=targetattrs, attrs2=outattrs, edgeLine=outattrs['lines']['line'])

        else:
            mgraph.addEdge(targetname, outname, attrs1=targetattrs, attrs2=outattrs, edgeLine=outattrs['lines']['line'])
            #mgraph.addEdge(outname, targetname, attrs1=outattrs, attrs2=targetattrs)

        return None

    elif functionTree.root.isParen:

        outname = output + '__' + location['subprogs'][-1]
        outattrs = {'location': location['subprogs'], 'cname': output, 'lines': location['lines']}
        #print targetname, outname

        if isOut:
            mgraph.addEdge(outname, targetname, attrs1=outattrs, attrs2=targetattrs, edgeLine=outattrs['lines']['line'])
            #mgraph.addEdge(targetname, outname, attrs1=targetattrs, attrs2=outattrs)

        elif 'direction' in targetattrs:
            #print targetattrs['direction']
            if targetattrs['direction'] == 'in':
                mgraph.addEdge(outname, targetname, attrs1=outattrs, attrs2=targetattrs, edgeLine=outattrs['lines']['line'])

            elif targetattrs['direction'] == 'out':
                mgraph.addEdge(targetname, outname, attrs1=targetattrs, attrs2=outattrs, edgeLine=outattrs['lines']['line'])

            # is 'inout'
            else:
                mgraph.addEdge(outname, targetname, attrs1=outattrs, attrs2=targetattrs, edgeLine=outattrs['lines']['line'])
                mgraph.addEdge(targetname, outname, attrs1=targetattrs, attrs2=outattrs, edgeLine=outattrs['lines']['line'])

        else:
            # mgraph.addEdge(targetname, outname, attrs1=targetattrs, attrs2=outattrs)
            mgraph.addEdge(outname, targetname, attrs1=outattrs, attrs2=targetattrs, edgeLine=outattrs['lines']['line'])

        return None

    else:

        raise RuntimeError("Unhandled case")



def rootVars(mgraph, instring, functions, location):
    # location = copy.deepcopy(location)
    stmtString = instring
    stmtString = ''.join(stmtString.split())
    sortfuncs = sorted(list(functions), key=len, reverse=True)

    for i in sortfuncs:
        tmp = ''.join(i.split())
        tmpname = tmp.split('(', 1)[0]
        stmtString = stmtString.replace(tmp, tmpname)

    rFuncs, _ = stmt_to_fnnames(stmtString)
    remFuncs = set([f.split('(', 1)[0] for f in functions])
    notFuncs = rFuncs - remFuncs

    rootvars = []
    for r in notFuncs:
        rootname = r + '__' + location['subprogs'][-1]
        rootattrs = {'location': location['subprogs'], 'cname': r, 'lines': location['lines']}
        rootvars.append((rootname, rootattrs))

    return rootvars


def deepStmts(mgraph, string, location, targetname, targetattrs, isOut=True):
    # location = copy.deepcopy(location)
    # targetattrs = copy.deepcopy(targetattrs)

    stringNoSpace = ''.join(string.split())
    _, functions = stmt_to_fnnames(stringNoSpace)

    # Need this complex parsing for statements like: alpha(beta) + foo(bar, alpha(beta))
    while functions:
        toRemove = []
        for i in functions:
            isRoot = True
            for j in functions:
                for k in parenthetic_contents(j):
                    substring = k[1]
                    if i in substring:
                        isRoot = False
                        break

            if isRoot:
                # It's a tree root, add edge to target.
                name = i.split('(', 1)[0]
                functionTree = FunctionTree(name, i, location, mgraph)
                functionTree.connectArgTree()
                connectOutput(functionTree, location, mgraph, targetname=targetname, targetattrs=targetattrs, isOut=isOut)

                # Remove from input string and continue loop.
                inoSpace = ''.join(i.split())
                stringNoSpace = stringNoSpace.replace(inoSpace, '', 1)
                toRemove.append(i)

        functions = functions.difference(toRemove)
        toRemove = []
        for f in functions:
            fnoSpace = ''.join(f.split())
            if fnoSpace not in stringNoSpace:
                toRemove.append(f)
        functions = functions.difference(toRemove)

    return None


def processSubArgs(mgraph, carg, defarg, location, subLoc):
    #carg = copy.deepcopy(carg)
    #defarg = copy.deepcopy(defarg)
    #location = copy.deepcopy(location)
    #subLoc = copy.deepcopy(subLoc)
    dtmp = ''.join(defarg[0].strip().split())
    defargName = dtmp + '__' + subLoc[-1]
    lines = mgraph.node[defargName]['attrs']['lines']
    defargAttrs = {'location': subLoc, 'cname': dtmp, 'direction': defarg[1], 'lines': lines}

    # subeq or function call
    if '=' in carg[0]:

        splits = carg[0].split('=', 1)
        if '%' not in splits[0]:
            lhs = splits[0].strip().lower().split('(')[0]
            lhs = ''.join(lhs.split())

        else:
            lhs = re.sub(r'\([^)]*\)', '', splits[0].strip().lower())
            # Sometimes this leaves a trailing ')'
            if ')' in lhs:
                lhs = lhs.replace(')', '')
                
            lhs = ''.join(lhs.split())

        lhsattrs = {'location': location['subprogs'], 'cname': lhs, 'lines': location['lines']}
        lhsName = lhs + '__' + location['subprogs'][-1]

        rhs = splits[1].strip().lower()

        names, functions = stmt_to_fnnames(rhs)
        rVars = rootVars(mgraph, rhs, functions, location)

        if defarg[1] == 'out':
            for rv in rVars:
                #mgraph.addEdge(rv[0], defargName, attrs1=rv[1], attrs2=defargAttrs)
                mgraph.addEdge(defargName, rv[0], attrs1=defargAttrs, attrs2=rv[1], edgeLine=rv[1]['lines']['line'])
                #print 'equal rv out', defargName, rv[0], defargAttrs, rv[1]

            #mgraph.addEdge(lhs, defargName, attrs1=defargAttrs, attrs2=defargAttrs)
            mgraph.addEdge(defargName, lhsName, attrs1=defargAttrs, attrs2=lhsattrs, edgeLine=lhsattrs['lines']['line'])

            #print 'out and equal', defargName, lhsName, defargAttrs, lhsattrs

            if functions:
                deepStmts(mgraph, rhs, location, defargName, defargAttrs, isOut=False)

        else:
            for rv in rVars:
                #mgraph.addEdge(defargName, rv[0], attrs1=defargAttrs, attrs2=rv[1])
                mgraph.addEdge(rv[0], lhsName, attrs1=rv[1], attrs2=lhsattrs, edgeLine=lhsattrs['lines']['line'])
                #print 'equal rv not out', rv[0], lhsName, rv[1], lhsattrs

            #mgraph.addEdge(defargName, lhs, attrs1=defargAttrs, attrs2=defargAttrs)
            mgraph.addEdge(lhsName, defargName, attrs1=lhsattrs, attrs2=defargAttrs, edgeLine=lhsattrs['lines']['line'])
            #print 'else and equal', lhsName, defargName, defargAttrs, lhsattrs

            if functions:
                deepStmts(mgraph, rhs, location, defargName, defargAttrs, isOut=False)

    else:

        names, functions = stmt_to_fnnames(carg[0])
        rVars = rootVars(mgraph, carg[0], functions, location)

        if defarg[1] == 'out':
            for rv in rVars:
                #mgraph.addEdge(rv[0], defargName, attrs1=rv[1], attrs2=defargAttrs)
                mgraph.addEdge(defargName, rv[0], attrs1=defargAttrs, attrs2=rv[1], edgeLine=rv[1]['lines']['line'])
                #print 'rv out', defargName, rv[0], defargAttrs, rv[1]

            if functions:
                deepStmts(mgraph, carg[0], location, defargName, defargAttrs, isOut=False)

        else:
            if defarg[1] == 'in':
                for rv in rVars:
                    #mgraph.addEdge(defargName, rv[0], attrs1=defargAttrs, attrs2=rv[1])
                    mgraph.addEdge(rv[0], defargName, attrs1=rv[1], attrs2=defargAttrs, edgeLine=rv[1]['lines']['line'])
                    #print 'rv in', rv[0], defargName, rv[1], defargAttrs

            elif defarg[1] == 'inout':
                for rv in rVars:
                    mgraph.addEdge(defargName, rv[0], attrs1=defargAttrs, attrs2=rv[1], edgeLine=rv[1]['lines']['line'])
                    mgraph.addEdge(rv[0], defargName, attrs1=rv[1], attrs2=defargAttrs, edgeLine=rv[1]['lines']['line'])
                    #print 'rv inout', rv[0], defargName, rv[1], defargAttrs                

            if functions:
                deepStmts(mgraph, carg[0], location, defargName, defargAttrs, isOut=False)

    return None
