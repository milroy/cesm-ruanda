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
from sklearn.model_selection import GridSearchCV
from argparse import ArgumentParser
import os.path

parser = ArgumentParser()
parser.add_argument("--test", dest="testfile", required=True, type=str)
parser.add_argument("--ensemble", dest="ensemblefile", required=True, type=str)
parser.add_argument("--ensSumUF", dest="ensemblefileUF", required=True, type=str)
parser.add_argument("--ensSumYR", dest="ensemblefileYR", required=True, type=str)
parser.add_argument("--timestep", dest="ts", required=True, type=int)
parser.add_argument("--nRuns", dest="num", required=True, type=int)
parser.add_argument("--maxRegCoef", dest="maxRegCoef", type=float)
parser.add_argument("--minRegCoef", dest="minRegCoef", type=float)
parser.add_argument("--stepSize", dest="stepSize", type=float)
parser.add_argument("--uf", dest="isUF", action="store_true")
parser.add_argument("--yr", dest="isYR", action="store_true")
parser.add_argument(
    "--standardize", dest="standard", action="store_true", default=False
)
parser.add_argument("--normalize", dest="normal", action="store_true", default=False)

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
maxRegCoef = args.maxRegCoef
minRegCoef = args.minRegCoef
stepSize = args.stepSize

regCoefs = np.arange(minRegCoef, maxRegCoef, stepSize)

print("Regularization coefficient values to test:")
print(regCoefs)

baseUF = nc.Dataset(ensemblefileUF, "r")

print(testfile, ts, num, isUF, isYR)

varssUF = []
for i in range(baseUF.variables["vars"].shape[0]):
    varssUF.append(baseUF.variables["vars"][i, :].tostring().strip())

baseYr = nc.Dataset(ensemblefileYR, "r")
varssYr = []
for i in range(baseYr.variables["vars"].shape[0]):
    varssYr.append(baseYr.variables["vars"][i, :].tostring().strip())


def build_examples(testSet, enssSet, ts, num, standard=False, normal=False):
    testArr = np.array(testSet, dtype=np.float64)[:num, :, ts]
    # testArr = testArr[np.random.choice(testArr.shape[0], num, replace=False), :]
    exSized = np.array(enssSet, dtype=np.float64)[:num, :, ts]
    # exSized = exSized[np.random.choice(exSized.shape[0], num, replace=False), :]
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
    labels = np.array(len(testArr) * [0.0] + exSized.shape[0] * [1.0])
    rng_state = np.random.get_state()
    np.random.shuffle(exs)
    np.random.set_state(rng_state)
    np.random.shuffle(labels)
    return exs, labels


with open(ensemblefile, "rb") as f:
    ensDict = pickle.load(f, encoding="latin1")

if not isinstance(ensDict, dict):
    ensDict = dict(ensDict)

const = [
    "BURDENSEASALT",
    "BURDENDUST",
    "BURDENSOA",
    "BURDENSO4",
    "BURDENPOM",
    "BURDENBC",
    "AODDUST3",
    "AODDUST1",
    "AODVIS",
]

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

with open(testfile, "rb") as f:
    testDict = pickle.load(f, encoding="latin1")

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

examples, labels = build_examples(
    testset, ensembleset, ts, num, standard=standard, normal=normal
)

print("input shape: ", examples.shape)
print("input rank: ", np.linalg.matrix_rank(examples))
print("condition: ", np.linalg.cond(examples))

parameters = {"C": regCoefs}

lr1 = linear_model.LogisticRegression(
    n_jobs=1, penalty="l1", solver="saga", max_iter=10000
)

clf = GridSearchCV(lr1, parameters, cv=5, n_jobs=3)
clf.fit(examples, labels)

print("Cross validation results:")
print(clf.cv_results_)

print("Cross validation best parameters:")
print(clf.best_params_)

print("Cross validation best estimator:")
print(clf.best_estimator_)
