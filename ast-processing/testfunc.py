def get_names(node, bag, f2003_type, depth):
    if isinstance(node, f2003_type):# and node.string not in bag:
        if isinstance(bag, set):
            bag.add(node.string)
        elif isinstance(bag, deque):
            bag.append(node.string)

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

def gen_argtree(fstring):
    def parenthetic_contents(string):
        """
        From stackoverflow question 4284991:
        Generate parenthesized contents in string as pairs (level, contents).
        """
        stopChars = set([' ', ',', '+', '-', '*', '/', '^'])
        stack = []
        for i, c in enumerate(string):
            if c == '(':
                stack.append(i)
            elif c == ')' and stack:
                func = ''
                start = stack.pop()
                if not stack == []:
                    func = re.search(r"[\w']+", string[stack[-1] + 1:start]).group()
                output = string[start + 1: i]
                if not output == '':
                    first = re.findall(r"[\w']+", output)[0]
                else:
                    first, output = '', ''
                yield (len(stack), first, func, output)
    argTree = {}
    for rets in parenthetic_contents(fstring):
        if not rets[0] in argTree:
            argTree[rets[0]] = {rets[1]: [{'expression': rets[3], 'parent': rets[2], \
                                            'output': rets[3], 'inFromBelow' : rets[3]}]}               
        elif not rets[1] in argTree[rets[0]]:
            argTree[rets[0]][rets[1]] = [{'expression': rets[3], 'parent': rets[2], \
                                                    'output': rets[3], 'inFromBelow' : rets[3]}]
        else:
            argTree[rets[0]][rets[1]].append({'expression': rets[3], 'parent': rets[2], \
                                                    'output': rets[3], 'inFromBelow' : rets[3]})
    return argTree


def processArgs(string):
    processedArgs = deque()
    string = string.replace(',', '+').replace(':', '+')
    if not string == '':
        reader = FortranStringReader('INVALID = ' + string)
        stmt = Fortran2003.Assignment_Stmt(reader)
        traverse(stmt, get_names, processedArgs, fparser.Fortran2003.Name)
        #print string, stmt.tostr()
        processedArgs.popleft()
    else:
        processedArgs = ['']
    return processedArgs 

def funcLiketoDict(string, argtree, fdict, isroot):
    isFunc = False
    isSubr = False
    isArray = False
    defArgs = None
    if not string in fdict:
        fdict[string] = {}
    if not 'argTree' in fdict[string]:
        fdict[string]['argTree'] = [argtree]
    else:
        fdict[string]['argTree'].append(argtree)
    try:
        funcDef = CESM.functions[string]
        defArgs = funcDef.args
        out = funcDef.output
        isFunc = True
    except KeyError:
        try:
            subroutineDef = CESM.subroutines[string]
            defArgs = subroutineDef.args
            out = subroutineDef.args
            isSubr = True
        except KeyError:
            out = string
            isArray = True
            pass
    fdict[string]['defArgs'] = defArgs
    fdict[string]['output'] = out
    fdict[string]['isRoot'] = isroot
    fdict[string]['isFunc'] = isFunc
    fdict[string]['isSubr'] = isSubr
    fdict[string]['isArray'] = isArray


def connectSubprogram(subp, callArgs):
    definitionArgs = subp['defArgs']
    print subp, definitionArgs
    prcArgs = processArgs(callArgs)
    outputs = subp['output']
    if not definitionArgs is None:
        for call, defArg in zip(prcArgs, definitionArgs):
            CESM.addEdge(defArg, call)
            if subp['isSubr']:
                # add intent in/out parsing
                continue
    else:
        # is array
        for call in prcArgs:
            CESM.addEdge(outputs, call)
    if subp['isSubr']:
        outputs = prcArgs
    return outputs



def setAboveInputs(tree, instring, outstring, depth):
    if depth > 0:
        for idx in xrange(len(tree[depth -1])):
            #print tree[depth -1], idx
            if tree[depth -1][idx]['stInput'].find(instring) > -1:
                #print idx, tree[depth -1][idx]
                tree[depth -1][idx]['inFromBelow'] = tree[depth -1][idx]['inFromBelow'].replace(instring, outstring)


def connectInOut(tree, depth):
    levelOutputs = []
    for idx, expr in enumerate(tree[depth]):
        #print idx, depth, expr
        searchExpr = expr['expression']
        if searchExpr in functions:
            procExpr = searchExpr.split('(', 1)
            #out = funcLike[procExpr[0]]['output']
            if funcLike[procExpr[0]]['isSubr']:
                # Connect edges and set output
                outProc = connectSubprogram(funcLike[procExpr[0]], procExpr[1][:-1])
                out = '+'.join(outProc)
                searchExpr = 'call(' + searchExpr + ')'
            else:
                out = connectSubprogram(funcLike[procExpr[0]], procExpr[1][:-1])
            expr['output'] = out
        elif expr['inFromBelow'] in functions:
            procExpr = expr['inFromBelow'].split('(', 1)
            #out = funcLike[procExpr[0]]['output']
            if funcLike[procExpr[0]]['isSubr']:
                # Connect edges and set output
                outProc = connectSubprogram(funcLike[procExpr[0]], procExpr[1][:-1])
                out = '+'.join(outProc)
                searchExpr = 'call(' + expr['inFromBelow'] + ')'
            #print depth, out
            else:
                out = connectSubprogram(funcLike[procExpr[0]], procExpr[1][:-1])
            expr['output'] = out
        else:
            for subexpr in expr['inFromBelow'].split(','):
                subexpr = subexpr.strip()
                #print subexpr
                if subexpr in functions:
                    procExpr = subexpr.split('(', 1)
                    #out = funcLike[procExpr[0]]['output']
                    if funcLike[procExpr[0]]['isSubr']:
                        # Connect edges and set output
                        outProc = connectSubprogram(funcLike[procExpr[0]], procExpr[1][:-1])
                        out = '+'.join(outProc)
                        searchExpr = 'call(' + subexpr + ')'
                    else:
                        out = connectSubprogram(funcLike[procExpr[0]], procExpr[1][:-1])
                    expr['inFromBelow'] = expr['inFromBelow'].replace(subexpr, out)
            expr['output'] = expr['inFromBelow']
        levelOutputs.append((searchExpr, expr['output']))
    return levelOutputs

def traverseTree(tree, depth):
    while depth > 0:
        levelOutputs = connectInOut(tree, depth)
        #print levelOutputs
        for i in levelOutputs:
            setAboveInputs(tree, i[0], i[1], depth)
        depth -= 1
    levelOutputs = connectInOut(tree, depth)
    for m in levelOutputs:
        setAboveInputs(tree, m[0], m[1], depth)

import fparser
from fparser import api
from fparser.api import Fortran2003
from fparser.readfortran import FortranStringReader
# Now the real work
s = 'rc = alpha(iota) * get_calday() * alpha(fubar) + delta ** delta + 1._r8 / ((1._r8 / rsmx(ispec)) + (1._r8 / rlux(ispec)) + (1._r8 / (rdc + rclx(ispec))) + (1._r8 / (rac(index_season, wesveg) + rgsx(ispec))) + alpha(beta(gamma, delta, call(urbanradiation(epsilon + (zeta(foxtrot + echo)), hit(assoc), nclex(omega : yotta), lambda(:, 2:1 + magma*eta(gfunc), :))))))'
reader = FortranStringReader(s)
stmt = Fortran2003.Assignment_Stmt(reader)
lhs = stmt.items[0].string

names = set()
#argLists = set()
emptyFunc = set()
functions = set()
arraysLike = set()
#dupParens = set()
traverse(stmt, get_names, names, fparser.Fortran2003.Name)
#traverse(stmt, get_names, argLists, fparser.Fortran2003.Section_Subscript_List)
traverse(stmt, get_names, emptyFunc, fparser.Fortran2003.Structure_Constructor)
traverse(stmt, get_names, functions, fparser.Fortran2003.Part_Ref)
#traverse(stmt, get_names, dupParens, fparser.Fortran2003.Parenthesis)
traverse(stmt, get_names, arraysLike, fparser.Fortran2003.Array_Section)

functions.update(emptyFunc)
functions.update(arraysLike)

# Connect level 0 variables to lhs
stmtString = stmt.tostr()
sfuncs = sorted(list(functions))
for i in sfuncs:
    stmtString = stmtString.replace(i, '')

for p in stmtString.split():
    if p in names and not p == lhs:
        CESM.addEdge(lhs, p)

# First make a mapping between functions and arrays and their args.
# Map and delete tree depth 0 functions
funcLike = {}
for i in functions:
    argTree = gen_argtree(i)
    isRoot = False
    depth = 0
    for j in functions:
        if i in j:
            depth += 1
    if depth == 1:
        isRoot = True
    name = i.split('(', 1)[0]
    funcLiketoDict(name, argTree, funcLike, isRoot)
    #print i, isSubr
    # It's a tree root, add edge to lhs.
    if isRoot and len(argTree) == 1:
        if not funcLike[name]['isArray']:
            definitionArgs = funcLike[name]['defArgs']
            # Know that argTree is depth 0, and the inputargs are the last added.
            inputArgs = funcLike[name]['argTree'][-1][0][-1]['expression']
            procArgs = processArgs(inputArgs)
            for call, defa in zip(procArgs, definitionArgs):
                for subCall in call:
                    CESM.addEdge(defa, subCall)
                    print defa, subCall
            if not funcLike[name]['isSubr']: 
                CESM.addEdge(lhs, funcLike[name]['output'])
                print lhs, funcLike[name]['output']
            else:
                for dArg in definitionArgs:
                    CESM.addEdge(lhs, dArg)
        # Must be array
        else:
            CESM.addEdge(lhs, funcLike[name]['output'])

            #print i, lhs, funcLike[name]['output']
    #elif isRoot: print name, funcLike[name]


for obj in funcLike:
    if funcLike[obj]['isRoot']:
        for tree in funcLike[obj]['argTree']:
            # root and deep
            depth = len(tree)
            if depth > 1:
                traverseTree(tree, depth -1)
