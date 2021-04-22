#!/bin/python

import sys
import os
import re
import copy
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
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--src", dest="srcdir", required=True, type=str)

args = parser.parse_args()
srcdir = args.srcdir

regex = r'\w+'
removeWords = ['.and.', '.or.']
likeFunc = re.compile("[\w']+\(+")

def getnames(stmt):

    names, functions = stmt_to_fnnames(stmt)

    hasFuncs = False
    if functions:
        hasFuncs = True

    return names, hasFuncs


def parse_eq(string, regex, removeW):
    # TODO: use F2003 parser
    tmp = re.findall(regex, string)
    resultWords  = [word for word in tmp if word not in removeW]

    return resultWords


def process_stmt(stmt, graphstack, nstack):
    namestack = copy.deepcopy(nstack)

    if isinstance(stmt, fparser.block_statements.Module):
        modname = stmt.name.lower()
        modVars = {''.join(x.translate(None, ''.join(["'", '/', '"'])).strip().lower().split()) 
                   for x in stmt.a.variable_names}
        modPublic = {''.join(x.translate(None, ''.join(["'", '/', '"'])).strip().lower().split()) 
                     for x in stmt.a.public_id_list}
        modPrivate = {''.join(x.translate(None, ''.join(["'", '/', '"'])).strip().lower().split()) 
                      for x in stmt.a.private_id_list}
        modProvides = {''.join(x.translate(None, ''.join(["'", '/', '"'])).strip().lower().split()) 
                       for x in stmt.a.module_provides.keys()}
        modSubprogs = {''.join(x.translate(None, ''.join(["'", '/', '"'])).strip().lower().split()) 
                       for x in stmt.a.module_subprogram.keys()}
        graphstack[0].modules[modname] = {'provides': modProvides, 'vars': modVars, 
                                          'public': modPublic, 'private': modPrivate, 
                                          'location': namestack, 'subprogs': modSubprogs}
        namestack.append(modname)

        lineNum = stmt.item.span[0]
        procNames = []
        for var in modVars:
            varattrs = {'location': namestack, 'cname': var, 'lines': {'line': lineNum}}
            varname = var + '__' + namestack[-1]
            procNames.append(varname)
            graphstack[0].addNode(varname, attrs=varattrs)

        graphstack.append(procNames)

    elif isinstance(stmt, fparser.block_statements.EndModule):
      # TODO: second parse with F2003 a la:
      # reader = FortranStringReader(stmt.tofortran())
      # subroutine = Fortran2003.End_Subroutine_Stmt(reader)
        subgraphNodes = graphstack.pop()
        module = graphstack[0].subgraph(subgraphNodes)
        graphstack[0].modules[stmt.name.lower()]['subgraph'] = module
        namestack.pop()

    elif isinstance(stmt, fparser.statements.Call):
      # TODO: second parse with F2003 a la:
      # reader = FortranStringReader(stmt.tofortran())
      # call = Fortran2003.Call_Stmt(reader)
        callName = stmt.designator.translate(None, ''.join(["'", '/', '"'])).lower()
        callArgs = [''.join(x.translate(None, ''.join(["'", '/', '"'])).strip().lower().split()) 
                    for x in stmt.items]
        lineNum = stmt.item.span[0]
        parentType = stmt.parent.blocktype
        parentLine = stmt.parent.item.span[0]
        lineInfo = {'line': lineNum, 'parentType': parentType, 'parentLine': parentLine}
        location = {'subprogs': namestack, 'lines': lineInfo}

        caList = []
        for v, ca in zip(stmt.items, callArgs):
            var_type = 'other'
            if stmt.get_variable(v).is_intent_in():
                var_type = 'in'
            elif stmt.get_variable(v).is_intent_out():
                var_type = 'out'
            elif stmt.get_variable(v).is_intent_inout():
                var_type = 'inout'

            caList.append((ca, var_type))

        if not callArgs == []:
            callArgs_0 = ''.join(callArgs[0].translate(None, ''.join(["'", '/', '"', '.', ',', '!', ':']))\
                       .strip().lower().split())

            # Adding 'processed': False for now, as connectivity will be affected.
            if callName == 'addfld':
                graphstack[0].calls.append((callName, caList, location, {'processed': False}))
              
                name1, _ = getnames(callArgs_0)

                for n1 in name1:
                    attrs1 = {'location': location['subprogs'], 'cname': n1, 'input': True, 'lines': location['lines']}
                    name1loc = n1 + '__' + location['subprogs'][-1]
                    graphstack[-1].append(name1loc)
                    graphstack[0].addNode(name1loc, attrs=attrs1)
                    graphstack[0].inputs.add(name1loc)

            elif callName == 'outfld':
                graphstack[0].calls.append((callName, caList, location, {'processed': False}))
                callArgs_1 = ''.join(callArgs[1].translate(None, ''.join(["'", '/', '"', '.', ',', '!', ':']))\
                           .strip().lower().split())

                name1, _ = getnames(callArgs_0)
                namen2, _ = getnames(callArgs_1)
                        
                for n1 in name1:
                    attrs1 = {'location': location['subprogs'], 'cname': n1, 'output': True, 'lines': location['lines']}
                    name1loc = n1 + '__' + location['subprogs'][-1]

                    for n2 in namen2:
                        namen2loc = n2 + '__' + location['subprogs'][-1]
                        attrsn2 = {'location': location['subprogs'], 'cname': n2, 'output': True, 'lines': location['lines']}
                        graphstack[0].addEdge(namen2loc, name1loc, attrs1=attrsn2, attrs2=attrs1, edgeLine=attrs1['lines']['line'])
                        graphstack[0].outputs.add(name1loc)
                        graphstack[0].outputs.add(namen2loc)
                        graphstack[-1].append(namen2loc)

            elif callName == 'infld':
                graphstack[0].calls.append((callName, caList, location, {'processed': False}))

                name1, _ = getnames(callArgs_0)

                if len(callArgs) == 15:
                    namen2, _ = getnames(callArgs[11])

                elif len(callArgs) == 14:
                    namen2, _ = getnames(callArgs[11])

                elif len(callArgs) == 12:
                    namen2, _ = getnames(callArgs[8])

                elif len(callArgs) == 11:
                    namen2, _ = getnames(callArgs[8])

                elif len(callArgs) == 7:
                    namen2, _ = getnames(callArgs[5])

                elif len(callArgs) == 5:
                    namen2, _ = getnames(callArgs[3])

                else:
                    print callArgs, namestack
                    print "Error processing infld call"
                    sys.exit(1)

                for n1 in name1:
                    attrs1 = {'location': location['subprogs'], 'cname': n1, 'input': True, 'lines': location['lines']}
                    name1loc = n1 + '__' + location['subprogs'][-1]

                    for n2 in namen2:
                        namen2loc = n2 + '__' + location['subprogs'][-1]
                        attrsn2 = {'location': location['subprogs'], 'cname': n2, 'input': True, 'lines': location['lines']}
                        graphstack[0].addEdge(namen2loc, name1loc, attrs1=attrsn2, attrs2=attrs1, edgeLine=attrs1['lines']['line'])
                        graphstack[0].inputs.add(name1loc)
                        graphstack[0].inputs.add(namen2loc)
                        graphstack[-1].append(namen2loc)

            else:
                graphstack[0].calls.append((callName, caList, location, {'processed': False}))

    elif isinstance(stmt, fparser.block_statements.Interface):
        ifacename = stmt.name.lower()
        procedures = [''.join(x.translate(None, ''.join(["'", '/', '"'])).strip().lower().split()) 
                      for x in stmt.a.module_procedures]

        if ifacename:
            graphstack[0].interfaces.append((ifacename, procedures))

    elif isinstance(stmt, fparser.block_statements.Subroutine):
      # TODO: second parse with F2003 a la:
      # reader = FortranStringReader(stmt.tofortran())
      # subroutine = Fortran2003.Subroutine_Stmt(reader)
        subname = stmt.name.lower()
        subArgs = [''.join(x.translate(None, ''.join(["'", '/', '"'])).strip().lower().split()) 
                   for x in stmt.args]
        saList = []
        for v, sa in zip(stmt.args, subArgs):
            var_type = 'other'
            if stmt.a.variables[v].is_intent_in():
                var_type = 'in'
            elif stmt.a.variables[v].is_intent_out():
                var_type = 'out'
            elif stmt.a.variables[v].is_intent_inout():
                var_type = 'inout'

            saList.append((sa, var_type))

        subVars = [''.join(x.translate(None, ''.join(["'", '/', '"'])).strip().lower().split()) 
                   for x in stmt.a.variables]
        graphstack[0].subroutines[subname] = {'args': saList, 'vars': subVars, 'location': namestack,
                                              'isIface': False}
        namestack.append(subname)

        lineNum = stmt.item.span[0]
        procNames = []
        for var in subVars:
            varattrs = {'location': namestack, 'cname': var, 'lines': {'line': lineNum}}
            varname = var + '__' + namestack[-1]
            graphstack[0].addNode(varname, attrs=varattrs)
            procNames.append(varname)

        graphstack.append(procNames)
        
    elif isinstance(stmt, fparser.block_statements.EndSubroutine):
      # TODO: second parse with F2003 a la:
      # reader = FortranStringReader(stmt.tofortran())
      # subroutine = Fortran2003.End_Subroutine_Stmt(reader)
        subgraphNodes = graphstack.pop()
        subroutine = graphstack[0].subgraph(subgraphNodes)
        graphstack[0].subroutines[stmt.name.lower()]['subgraph'] = subroutine
        namestack.pop()

    elif isinstance(stmt, fparser.block_statements.Function):
      # TODO: second parse with F2003 a la:
      # reader = FortranStringReader(stmt.tofortran())
      # function = Fortran2003.Function_Stmt(reader)
        fname = stmt.name.lower()
        funcArgs = [''.join(x.translate(None, ''.join(["'", '/', '"'])).strip().lower().split()) 
                    for x in stmt.args]
        output = stmt.result
        funcVars = [''.join(x.translate(None, ''.join(["'", '/', '"'])).strip().lower().split()) 
                    for x in stmt.a.variables]
        graphstack[0].functions[fname] = {'args': funcArgs, 'vars': funcVars, 'location': namestack, 
                                          'output': output, 'isIface': False, 'isIntrinsic': False}
        namestack.append(fname)

        lineNum = stmt.item.span[0]
        procNames = []
        for var in funcVars:
            varattrs = {'location': namestack, 'cname': var, 'lines': {'line': lineNum}}
            varname = var + '__' + namestack[-1]
            graphstack[0].addNode(varname, attrs=varattrs)
            procNames.append(varname)

        graphstack.append(procNames)

    elif isinstance(stmt, fparser.block_statements.EndFunction):
      # TODO: second parse with F2003 a la:
      # reader = FortranStringReader(stmt.tofortran())
      # function = Fortran2003.End_Function_Stmt(reader)
        subgraphNodes = graphstack.pop()
        function = graphstack[0].subgraph(subgraphNodes)
        graphstack[0].functions[stmt.name.lower()]['subgraph'] = function
        namestack.pop()

    elif isinstance(stmt, fparser.block_statements.Use):
      # TODO: second parse with F2003 a la:
      # reader = FortranStringReader(stmt.tofortran())
      # use = Fortran2003.Use_Stmt(reader)
        renames = {'same': [], 'renamed': [], 'all': False}

        if stmt.items:
            for item in stmt.items:

                splits = item.strip().lower().replace(' ', '').split('=>')
                # Fortran use rename is confusing: newName => sourceName
                if len(splits) < 2:
                    renames['same'].append(splits[0])

                else:
                    renames['renamed'].append((splits[0], splits[1]))

        else:
            renames['all'] = True

        graphstack[0].uses.append((stmt.name.lower(), renames, namestack))

    elif isinstance(stmt, fparser.block_statements.Assignment):
        if '%' not in stmt.variable:
            assigned = ''.join(stmt.variable.strip().lower().split()).split('(', 1)[0]
            lhs = assigned + '__' + namestack[-1]

        else:
            assigned = re.sub(r'\([^)]*\)', '', ''.join(stmt.variable.strip().lower().split()))
            if ')' in assigned:
                assigned = assigned.replace(')', '')

            lhs = assigned + '__' + namestack[-1]
            
        lineNum = stmt.item.span[0]
        parentType = stmt.parent.blocktype
        parentLine = stmt.parent.item.span[0]
        lineInfo = {'line': lineNum, 'parentType': parentType, 'parentLine': parentLine}
        location = {'subprogs': namestack, 'lines': lineInfo}
        aattrs = {'location': location['subprogs'], 'cname': assigned, 'lines': location['lines']}

        try:
            reader = FortranStringReader(stmt.tofortran())
            assnStmt = Fortran2003.Assignment_Stmt(reader)
            # Do some string processing
            rhs = assnStmt.items[2].tostr().strip().lower()
            # TODO: replace with F2003 parser
            #if likeFunc.search(rhs):
            names, hasFunc = getnames(rhs)

            if hasFunc:
                graphstack[0].funcLike.append((assnStmt, location))

            else:
                graphstack[-1].append(lhs)

                for j in names: 
                    jname = j + '__' + location['subprogs'][-1]
                    jattrs = {'location': location['subprogs'], 'cname': j, 'lines': location['lines']}
                    graphstack[0].addEdge(jname, lhs, attrs1=jattrs, attrs2=aattrs, edgeLine=aattrs['lines']['line'])
                    graphstack[-1].append(jname)
        except:
            try:
                reader = FortranStringReader(assigned + '=' + rhs)
                assnStmt = Fortran2003.Assignment_Stmt(reader)
                # Do some string processing
                rhs = assnStmt.items[2].tostr().strip().lower()
                names, hasFunc = getnames(rhs)
                if hasFunc:
                    graphstack[0].funcLike.append((assnStmt, location))

                else:
                    graphstack[-1].append(lhs)
                    for j in names: 
                        jname = j + '__' + location['subprogs'][-1]
                        jattrs = {'location': location['subprogs'], 'cname': j, 'lines': location['lines']}
                        graphstack[0].addEdge(jname, lhs, attrs1=jattrs, attrs2=aattrs, edgeLine=aattrs['lines']['line'])
                        graphstack[-1].append(jname)

            except:
                # Do some string processing
                rhs = ''.join(stmt.expr.strip().lower().split())
                rhs = re.sub(' +',' ', rhs)
                rhs = re.sub(' \(', '(', rhs)
                if likeFunc.search(rhs):
                    graphstack[0].funcLike.append((assigned + '=' + rhs, location))
                else:
                    graphstack[-1].append(lhs)
                    for j in parse_eq(rhs, regex, removeWords):
                        jname = j + '__' + location['subprogs'][-1]
                        jattrs = {'location': location['subprogs'], 'cname': j, 'lines': location['lines']}
                        graphstack[0].addEdge(jname, lhs, attrs1=jattrs, attrs2=aattrs, edgeLine=aattrs['lines']['line'])
                        graphstack[-1].append(jname)


    elif isinstance(stmt, fparser.block_statements.PointerAssignment):
        assignedSpace = stmt.variable.strip().lower()
        if '%' not in assignedSpace:
            assigned = ''.join(assignedSpace.split()).split('(', 1)[0]

        else:
            assigned = re.sub(r'\([^)]*\)', '', ''.join(assignedSpace.split()))
            if ')' in assigned:
                assigned = assigned.replace(')', '')

        lhs = assigned + '__' + namestack[-1]
        lineNum = stmt.item.span[0]
        parentType = stmt.parent.blocktype
        parentLine = stmt.parent.item.span[0]
        lineInfo = {'line': lineNum, 'parentType': parentType, 'parentLine': parentLine}
        location = {'subprogs': namestack, 'lines': lineInfo}
        aattrs = {'location': location['subprogs'], 'cname': assigned, 'lines': location['lines']}
        graphstack[-1].append(lhs)

        assigneeSpace = stmt.expr.strip().lower()
        if '%' not in assigneeSpace:
            assignee = ''.join(assigneeSpace.split()).split('(', 1)[0]

        else:
            assignee = re.sub(r'\([^)]*\)', '', ''.join(assigneeSpace.split()))
            if ')' in assignee:
                assignee = assignee.replace(')', '')

        rhs = assignee + '__' + location['subprogs'][-1]
        aaeettrs = {'location': location['subprogs'], 'cname': assignee, 'lines': location['lines']}
        graphstack[-1].append(rhs)

        graphstack[0].addEdge(lhs, rhs, attrs1=aattrs, attrs2=aaeettrs, edgeLine=aattrs['lines']['line'])
        graphstack[0].addEdge(rhs, lhs, attrs1=aaeettrs, attrs2=aattrs, edgeLine=aattrs['lines']['line'])

    return namestack

def readF(path):
    with open(path, 'r') as f:
        for line in f:
            yield line



CESM = metagraph('CESM')
stack = [CESM]
namestack = []
badfiles = []
for root, dirs, files in os.walk(srcdir):

    useFiles = [root + f for f in files]
    for file in files:
        if file.endswith('.F90') or file.endswith('.f90') or file.endswith('f.h'):
            try:
                namestack.append(file)
                tree = api.parse(os.path.join(root, file))#, source_only=useFiles)

                for stmt, depth in api.walk(tree):
                    namestack = process_stmt(stmt, stack, namestack)                        

            except Exception as ex:
                badfiles.append((file, ex))
                pass

            namestack.pop()



for i in Intrinsic_Procedures:
   if i not in CESM.functions:# and 'max' not in i and 'min' not in i:
       CESM.functions[i] = {'args': None, 'vars': None, 'location': ['CESM'], 
                            'output': i, 'isIface': False, 'isIntrinsic': True}


CESM.processIfaces()
CESM.processUsemaps()
ca = CESM.connectCallArgs()
cf = CESM.connectFunctions()
CESM.processUsedVars()
#a, b, c, d, e, eigen, comms = CESM.eigenCentrality(['freqr', 'freqs', 'cld', 'clmed', 'clhgh'], models=['cam'], hashi=True, girvan=True)