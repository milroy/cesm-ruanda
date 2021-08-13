#!/bin/python

import numpy as np
import glob
from collections import OrderedDict
from collections import defaultdict
import pickle
import itertools
import netCDF4 as nc
from sklearn.preprocessing import StandardScaler, Normalizer
from sklearn import linear_model
from argparse import ArgumentParser
import os.path

parser = ArgumentParser()
parser.add_argument("--test", dest="testfile", required=True, type=str)
parser.add_argument("--ensemble", dest="ensemblefile", required=True, type=str)
parser.add_argument("--ensSumUF", dest="ensemblefileUF", required=True, type=str)
parser.add_argument("--ensSumYR", dest="ensemblefileYR", required=True, type=str)
parser.add_argument("--timestep", dest="ts", required=True, type=int)
parser.add_argument("--nRuns", dest="num", required=True, type=int)
parser.add_argument("--regCoef", dest="regCoef", type=float)
parser.add_argument("--uf", dest="isUF", action='store_true')
parser.add_argument("--yr", dest="isYR", action='store_true')
parser.add_argument("--standardize", dest="standard", action='store_true', default=False)
parser.add_argument("--normalize", dest="normal", action='store_true', default=False)

args = parser.parse_args()
testfile = args.testfile
ensemblefile = args.ensemblefile
ensemblefileUF = args.ensemblefileUF
ensemblefileYR = args.ensemblefileYR
ts = args.ts
num = args.num
isUF = args.isUF
isYR = args.isYR
standard = args.standard
normal = args.normal
regCoef = args.regCoef

baseUF = nc.open_file(ensemblefileUF, 'r')

print(testfile, ts, num, isUF, isYR)

varssUF = []
for i in range(baseUF.variables['vars'].shape[0]):
    varssUF.append(baseUF.variables['vars'][i,:].tostring().strip())
    
baseYr = nc.open_file(ensemblefileYR, 'r')
varssYr = []
for i in range(baseYr.variables['vars'].shape[0]):
    varssYr.append(baseYr.variables['vars'][i,:].tostring().strip())

def build_examples(testSet, enssSet, ts, num, standard=False, normal=False):
    testArr = np.array(testSet, dtype=np.float64)[:num, :, ts]
    exSized = np.array(enssSet, dtype=np.float64)[:num, :, ts]
    if standard:
        ss = StandardScaler()
        exs = ss.fit_transform(np.vstack((testArr, exSized)))
    elif normal:
        ss = StandardScaler(with_std=False)
        nz = Normalizer()
        exs = ss.fit_transform(np.vstack((testArr, exSized)))
        exs = nz.fit_transform(exs.T).T
    else:
        exs = np.vstack((testArr, exSized))
    labels = np.array(len(testArr)*[0.] + exSized.shape[0]*[1.])
    rng_state = np.random.get_state()
    np.random.shuffle(exs)
    np.random.set_state(rng_state)
    np.random.shuffle(labels)
    return exs, labels


with open(ensemblefile, 'rb') as f:
    u = pickle._Unpickler(f)
    u.encoding = 'latin1'
    ensDict = u.load()

if not isinstance(ensDict, dict):
    ensDict = dict(ensDict)

const = ['BURDENSEASALT', 'BURDENDUST', 'BURDENSOA', 'BURDENSO4', 'BURDENPOM', 
         'BURDENBC', 'AODDUST3', 'AODDUST1', 'AODVIS']
    
ensemble = []
ensembleUF = []
for i in sorted(ensDict.keys()):
    ensemble.append(ensDict[i])
    tmp = []
    for j in range(ensDict[i].shape[0]):
        if varssYr[j] in varssUF:
            if not varssYr[j] in const:
                tmp.append(ensDict[i][j])
    ensembleUF.append(np.asarray(tmp))
    
with open(testfile, 'rb') as f:
    u = pickle._Unpickler(f)
    u.encoding = 'latin1'
    testDict = u.load()

if not isinstance(testDict, dict):
	testDict = dict(testDict)    

test = []
testUF = []
for i in sorted(testDict.keys()):
    test.append(testDict[i])
    tmp = []
    for j in range(testDict[i].shape[0]):
        if varssYr[j] in varssUF:
            if not varssYr[j] in const:
                tmp.append(testDict[i][j])
    testUF.append(np.asarray(tmp))

if isUF:
	testset = testUF
	ensembleset = ensembleUF
	Vs = varssUF

elif isYR:
	testset = test
	ensembleset = ensemble
	Vs = varssYr

examples, labels = build_examples(testset, ensembleset, ts, num, standard=standard, normal=normal)

print("input shape: ", examples.shape)
print("input rank: ", np.linalg.matrix_rank(examples))
print("condition number: ", np.linalg.cond(examples))

lr1 = linear_model.LogisticRegression(C=regCoef, n_jobs=4, penalty='l1', solver='saga', max_iter=10000)
lr1.fit(examples, labels)
betaslr1 = lr1.coef_[0]
sortlr1 = np.argsort(-np.abs(betaslr1))
svlr1p = [Vs[i].lower() for i in sortlr1 if betaslr1[i] > 0.]
svlr1n = [Vs[i].lower() for i in sortlr1 if betaslr1[i] < 0.]

nrawcoefs = np.array([betaslr1[i] for i in sortlr1 if betaslr1[i] < 0.])
prawcoefs = np.array([betaslr1[i] for i in sortlr1 if betaslr1[i] > 0.])

ensbs = zip(svlr1p, prawcoefs)
expbs = zip(svlr1n, nrawcoefs)

print("Ensemble betas: ", list(ensbs))
print("Experiment betas: ", list(expbs))
