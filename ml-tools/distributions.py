#!/usr/bin/python

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

def get_gms(base, test, path):
    vars_s = []
    for i in xrange(base.variables['vars'].shape[0]):
        vars_s.append(base.variables['vars'][i,:].tostring().strip())

    ensGms = base.variables['global_mean'][:].T.astype(np.float64)
    testGms = np.asarray(test).astype(np.float64)

    global_means = OrderedDict()
    global_means_sum = OrderedDict()
    files = sorted(glob.glob(path + '*'))
    for filename in files:
        filenumber = re.findall(r"\D(\d{3})\D", filename)[0]
        fn = int(filenumber)
        global_means_run = OrderedDict()
        global_means_run_sum = OrderedDict()

        i = 0
        for vari in vars_s:
            global_means_run_sum[vari] = ensGms[fn, i]
            global_means_run[vari] = testGms[fn, i]
            i += 1 

        global_means[fn] = global_means_run
        global_means_sum[fn] = global_means_run_sum


    return global_means, global_means_sum


def scale(ens_global_means, exp_global_means, nruns, tstep):
    diff_scaled_exp = OrderedDict()
    for sim in ens_global_means:
        if sim != 0 and sim < nruns:
            diffvar = OrderedDict()
            for var in ens_global_means[sim]:
                diffvar[var] = (exp_global_means[sim][var][tstep] - ens_global_means[0][var][tstep])/(ens_global_means[0][var][tstep])
                diff_scaled_exp[sim] = diffvar

    exp = defaultdict(OrderedDict)
    for sim in diff_scaled_exp:
        for var in diff_scaled_exp[sim]:
            exp[var][sim] = diff_scaled_exp[sim][var]

    diff_scaled_ens = OrderedDict()
    for sim in ens_global_means:
        if sim != 0 and sim < nruns:
            diffvar = OrderedDict()
            for var in ens_global_means[sim]:
                diffvar[var] = (ens_global_means[sim][var][tstep] - ens_global_means[0][var][tstep])/(ens_global_means[0][var][tstep])
            diff_scaled_ens[sim] = diffvar

    ens = defaultdict(OrderedDict)
    for sim in diff_scaled_ens:
        for var in diff_scaled_ens[sim]:
            ens[var][sim] = diff_scaled_ens[sim][var]

    return (ens, exp)


def classify(ens, exp, enslabel, explabel, tstep):

    case1 = []

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
            exp_l[:] = 0.

        if False in np.isfinite(ens_l):
            ens_l[:] = 0.

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
            case1.append((vary, diff))
        elif (exp_median > iqr_ens[0] or exp_median < iqr_ens[1]) and (ens_median > iqr_exp[0] or ens_median < iqr_exp[1]):
            pass
        else:
            pass
    return case1

def main(sumFile, ensFile, testFile, tstep):
    sumFile = str(sumFile)
    testFile = str(testFile)
    ensFile = str(ensFile)

    base = nc.Dataset(sumFile, 'r')
    vs = []
    for i in xrange(base.variables['vars'].shape[0]):
        vs.append(base.variables['vars'][i,:].tostring().strip())

    with open(ensFile, 'rb') as f:
        ens = pickle.load(f)

    ens_gms = defaultdict(OrderedDict)
    for sim in ens:
        for idx, var in enumerate(vs):
            ens_gms[sim][var] = ens[sim][idx]

    ens_gms = {int(key):ens_gms[key] for key in ens_gms}

    with open(testFile, 'rb') as f:
        testmeans = pickle.load(f)

    if isinstance(testmeans, list):
        testmeans = dict(testmeans)

    exp_gms = defaultdict(OrderedDict)
    for sim in testmeans:
        for idx, var in enumerate(vs):
            exp_gms[sim][var] = testmeans[sim][idx]


    exp_gms = {int(key):exp_gms[key] for key in exp_gms}

    numtestFiles = len(exp_gms)

    scaled_ens, scaled_exp = scale(ens_gms, exp_gms, numtestFiles, tstep)

    class1 = classify(ens_gms, exp_gms, 'ensemble', 'experiment', tstep)

    class1.sort(key=itemgetter(1), reverse=True)
    return class1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Dist script')
    parser.add_argument('sumFile', type=str, default='')
    parser.add_argument('testFile', type=str, default='')
    parser.add_argument('tstep', type=int, default=0)
    parser.add_argument('ensFile', type=str, default='')

    args = parser.parse_args()

    ordered = main(args.sumFile, args.ensFile, args.testFile, args.tstep)

