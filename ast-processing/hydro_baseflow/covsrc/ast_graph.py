import sys, os, re
sys.path.insert(0, '/home/milroy/git/f2py/')
import fparser
from fparser import api
from fparser.api import get_reader
from fparser.Fortran2003 import *
import networkx as nx
from collections import defaultdict, OrderedDict
from subgraphs import metagraph
#sys.path.insert(0, '/home/milroy/git/KGen/base/')
#import logging
##from numpy.distutils.misc_util import yellow_text, red_text # KGEN deletion
#from utils import split_comma, specs_split_comma, is_int_literal_constant
#from utils import classes
#import base_classes

rootdir = '/home/milroy/git/phd-code/scripts/cam-var-search/ast/src/'

removeChars = [')', '(', ',', ':', '+', '/', '-', '*', '==', '<', '>']
removeWords = ['.and.', '.or.']
rx = '[' + re.escape(''.join(removeChars)) + ']'

# def process_call(stmt, depth, callList):
#     var_type = {}
#     for v in stmt.items:
#         if stmt.get_variable(v).is_intent_in():
#             var_type[v] = 'in'
#         elif stmt.get_variable(v).is_intent_out():
#             var_type[v] = 'out'
#         elif stmt.get_variable(v).is_intent_inout():
#             var_type[v] = 'inout'
#     callList.append({'args': var_type, 'target': stmt.designator, 'depth': depth})
#     return None

# def process_subroutine(stmt, depth, subrDict):
#     subrDict[stmt.name] = {'args': stmt.args, 'vars': stmt.a.variables, 'depth': depth}
#     return None

# def process_function(stmt, depth, funcDict):
#     funcDict[stmt.name] = {'args': stmt.args, 'vars': stmt.a.variables, 'output': stmt.result, 'depth': depth}
#     return None

# def process_use(stmt, depth, useDict):
#     useDict[stmt.name] = {'only': stmt.items, 'depth': depth}
#     return None

def parse_eq(strng, removeX, removeW):
  tmp = re.sub(r'\([^)]*\)', '', strng.lower())
  tmp1 = re.sub(removeX, ' ', tmp).strip().split()
  resultWords  = [word for word in tmp1 if word not in removeW]
  #result = ' '.join(resultWords)
  #resultList = result.strip().split()
  # splits = strng.split('(')
  # for paren in splits:
  #   if ')' in paren:
  #       args = paren.split(')')[0]
  #       tmpobj = splits[splits.index(paren) - 1]
  #       obj = re.sub(removeX, ' ', tmpobj).strip().split()[-1]
  #       print obj, ' ;', type(args)
  #       refsDict[obj][assigned] = args
  return resultWords


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


def process_stmt(stmt, stack):
    if isinstance(stmt, fparser.statements.Call):
        stack[-1].calls[stmt.designator] = stmt.items

    elif isinstance(stmt, fparser.block_statements.Subroutine):
        subroutine = metagraph(stmt.name)
        subroutine.args = stmt.args
        for var in stmt.a.variables:
            subroutine.addNode(var)
        stack.append(subroutine)

    elif isinstance(stmt, fparser.block_statements.EndSubroutine):
        eSbr = stack.pop()
        stack[-1].addNode(eSbr)

    elif isinstance(stmt, fparser.block_statements.Function):
        function = metagraph(stmt.name)
        function.args = stmt.args
        for var in stmt.a.variables:
            function.addNode(var)
        stack.append(function)

    elif isinstance(stmt, fparser.block_statements.EndFunction):
        eFn = stack.pop()
        stack[-1].addNode(eFn)

    elif isinstance(stmt, fparser.block_statements.Use):
        stack[-1].uses[stmt.name] = stmt.items

    elif isinstance(stmt, fparser.block_statements.Assignment):
#        try:
#            for i in stmt.entity_decls:
#                tmp = i.strip().split()[0].lower().split('(', 1)[0]
#                stack[-1].addNode(tmp)
#        except:
#            pass
#        try:
        assigned = stmt.variable.lower().split('(', 1)[0]
        for j in parse_eq(stmt.expr, rx, removeWords):
            stack[-1].addEdge(j, assigned)



CESM = metagraph('CESM')
stack = [CESM]
for root, dirs, files in os.walk(rootdir):
     for file in files:
        if file.endswith('.F90_kgenout'):
            try:
                try:
                    try:
                        moduleGraph = metagraph(file)
                        tree = api.parse(os.path.join(root, file), source_only=files)
                        for stmt, depth in api.walk(tree):
                            #print stmt
                            process_stmt(stmt, stack)
                        stack[-1].addNode(moduleGraph)
                        stack.append(moduleGraph)
                        stack.pop()
                    except:
                        moduleGraph = metagraph(file)
                        tree = api.parse(os.path.join(root, file), source_only=None)
                        for stmt, depth in api.walk(tree):
                            #print stmt
                            process_stmt(stmt, stack)
                        stack[-1].addNode(moduleGraph)
                        stack.append(moduleGraph)
                        stack.pop()                        
                except:
                    moduleGraph = metagraph(file)
                    reader = get_reader(os.path.join(root, file))
                    prog = Program(reader)
                    tree = api.parse(prog)
                    for stmt, depth in api.walk(tree):
                        #print stmt
                        process_stmt(stmt, stack)
                    stack[-1].addNode(moduleGraph)
                    stack.append(moduleGraph)
                    stack.pop()

            except:
                print file
                pass




# for call in subProgs['call']:
#     sbr = call['target']
#     for arg in call['args']:
#         src = arg.lower()
#         idx = call['args'].index(arg)
#         try:
#             trgt = subProgs['subroutine'][sbr]['args'][idx].lower().split('(', 1)[0]
#             G.add_nodes_from([src, trgt])
#             G.add_edge(src, trgt)
#         except:
#             pass

# for func in subProgs['function']:
#     for var in refs[func]:
#         for call, arg in zip(refs[func][var], subProgs['function'][func]['args']):
#             G.add_edge(call, arg)
#         G.add_edge(subProgs['function'][func]['output'], var)
