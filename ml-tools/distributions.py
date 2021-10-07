#!/usr/bin/python3

import numpy as np
from collections import OrderedDict
from collections import defaultdict
import pickle
import itertools
import netCDF4 as nc
import glob
import re
import os
import argparse
from operator import itemgetter
from sklearn.preprocessing import StandardScaler


def classify(ens, exp, enslabel, explabel, tstep):

    class1 = []

    ensd = defaultdict(list)
    expd = defaultdict(list)
    for sim in exp:
        if sim in ens:
            for vary in ens[sim]:
                ensd[vary].append(ens[sim][vary][tstep])
                expd[vary].append(exp[sim][vary][tstep])

    for vary in ensd:
        exp_l = np.array(expd[vary]).reshape(-1, 1)
        ens_l = np.array(ensd[vary]).reshape(-1, 1)

        if False in np.isfinite(exp_l):
            exp_l[:] = 0.0

        if False in np.isfinite(ens_l):
            ens_l[:] = 0.0

        ss = StandardScaler()
        ss.fit(ens_l)
        ens_l = ss.transform(ens_l)
        exp_l = ss.transform(exp_l)

        exp_median = np.median(exp_l)
        ens_median = np.median(ens_l)
        iqr_exp = np.percentile(exp_l, [75, 25])
        iqr_ens = np.percentile(ens_l, [75, 25])

        if iqr_exp[1] > iqr_ens[0] or iqr_ens[1] > iqr_exp[0]:
            diff = abs(ens_median - exp_median)
            class1.append((vary, diff))

    return class1


def main(sumFile, ensFile, testFile, tstep):
    sumFile = str(sumFile)
    testFile = str(testFile)
    ensFile = str(ensFile)

    base = nc.Dataset(sumFile, "r")
    vs = []
    for i in range(base.variables["vars"].shape[0]):
        vs.append(base.variables["vars"][i, :].tostring().strip())

    with open(ensFile, "rb") as f:
        u = pickle._Unpickler(f)
        u.encoding = "latin1"
        ens = u.load()

    ens_gms = defaultdict(OrderedDict)
    for sim in ens:
        for idx, var in enumerate(vs):
            ens_gms[sim][var] = ens[sim][idx]

    ens_gms = {int(key): ens_gms[key] for key in ens_gms}

    with open(testFile, "rb") as f:
        u = pickle._Unpickler(f)
        u.encoding = "latin1"
        testmeans = u.load()

    if isinstance(testmeans, list):
        testmeans = dict(testmeans)

    exp_gms = defaultdict(OrderedDict)
    for sim in testmeans:
        for idx, var in enumerate(vs):
            exp_gms[sim][var] = testmeans[sim][idx]

    exp_gms = {int(key): exp_gms[key] for key in exp_gms}

    class1 = classify(ens_gms, exp_gms, "ensemble", "experiment", tstep)

    class1.sort(key=itemgetter(1), reverse=True)
    return class1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Dist script")
    parser.add_argument("sumFile", type=str, default="")
    parser.add_argument("testFile", type=str, default="")
    parser.add_argument("tstep", type=int, default=0)
    parser.add_argument("ensFile", type=str, default="")

    args = parser.parse_args()

    print(main(args.sumFile, args.ensFile, args.testFile, args.tstep))
