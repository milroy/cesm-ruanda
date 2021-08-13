#!/bin/python

import numpy as np
import glob
import os
import cPickle as pickle
import itertools
import netCDF4 as nc
from pyspark.sql import SparkSession
spark = SparkSession.builder.getOrCreate()
sc = spark.sparkContext
import pyspark
from pyspark import SparkConf, SparkContext
from sklearn import linear_model
import distributions, rlr_vs_bp
import numpy as np
from collections import OrderedDict, defaultdict
import seaborn as sns
from matplotlib.colors import ListedColormap
from sklearn.preprocessing import StandardScaler

partitions = 30
testDirYr = "/glade/p/tdd/asap/verification/cesm1_3_beta11/"
testDirUF = "/glade/p/tdd/asap/verification/UF_CECT/cesm1_3_beta11/"
testsYr = [""] 
testsUF = [""] 
testsUF = os.listdir(testDirUF)
baseYr = nc.Dataset('/glade/p/tdd/asap/verification/pca_ens_summary/cesm1_3_beta11/sz453.compilers-rand2_V6.nc', 'r')
baseUF = nc.Dataset('/glade/p/work/milroy/sz750-ufect_T.nc', 'r')

varssUF = []
for i in xrange(baseUF.variables['vars'].shape[0]):
    varssUF.append(baseUF.variables['vars'][i,:].tostring().strip())

varUF = sc.broadcast(varssUF)
gmsUF = baseUF.variables['global_mean'][:, :30].T.astype(np.float64)

varssYr = []
for i in xrange(baseYr.variables['vars'].shape[0]):
    varssYr.append(baseYr.variables['vars'][i,:].tostring().strip())

varYr = sc.broadcast(varssYr)
gmsYr = baseYr.variables['global_mean'][:, :30].T.astype(np.float64)

def gmsNetCDF(part):
    f = nc.Dataset(part, 'r')
    area = f.variables["area"]
    area_wgt = area[:].astype(np.float64)
    total = np.sum(area_wgt)
    area_wgt /= total
    nlev = f.dimensions["lev"]
    global_means_run = np.zeros((len(varUF.value), f.variables[varUF.value[0]].shape[0]), dtype=np.float64)
    k = 0
    for vari in varUF.value:
        if len(f.variables[vari].shape) == 3:
            per_tstep = np.zeros(f.variables[vari].shape[0], dtype=np.float64)
            for tstep in xrange(f.variables[vari].shape[0]):
                gm_lev = np.zeros(nlev, dtype=np.float64)
                for i in xrange(nlev):
                    gm_lev[i] = np.average(f.variables[vari][tstep, i, :].astype(np.float64), weights=area_wgt)
                    per_tstep[tstep] = np.mean(gm_lev)
            global_means_run[k] = per_tstep
        else:
            per_tstep = np.zeros(f.variables[vari].shape[0], dtype=np.float64)
            for tstep in xrange(f.variables[vari].shape[0]):
                per_tstep[tstep] = np.average(f.variables[vari][tstep, :].astype(np.float64), weights=area_wgt)
            global_means_run[k] = per_tstep
        k += 1
    return global_means_run


def build_examples(testSet, enssSet, nPerms, aug=True):
    if aug:
        ensSet = list(enssSet[:len(testSet)])
        exs = []
        for k in xrange(nPerms):
            np.random.shuffle(testSet)
            tmp = []
            for i in xrange((len(testSet))):
                tmp.append(testSet[i][:,-1])
            exs.append(np.array(tmp, dtype=np.float64).flatten())
        for k in xrange(nPerms):
            np.random.shuffle(ensSet)
            tmp = []
            for i in xrange(len(ensSet)):
                tmp.append(ensSet[i][:])
            exs.append(np.array(tmp, dtype=np.float64).flatten())
        labels = nPerms*[0.] + nPerms*[1.]
        labels = np.array(labels)
        exs = np.array(exs, dtype=np.float64)
        rng_state = np.random.get_state()
        np.random.shuffle(exs)
        np.random.set_state(rng_state)
        np.random.shuffle(labels)
    else:
        exs = np.vstack((np.array(testSet, dtype=np.float64)[:, :, -1], np.array(enssSet, dtype=np.float64)))
        labels = np.array(len(testSet)*[0.] + len(enssSet)*[1.])
        rng_state = np.random.get_state()
        np.random.shuffle(exs)
        np.random.set_state(rng_state)
        np.random.shuffle(labels)
    return exs, labels


#finished = ["clm_urban", "solar", "clm_default", "NU-P", "RH-MIN-LOW", "RAND_KV_DIFFSEED", "RAND-MT", "CONV_WATER"]
#finished = ["PREC", "RAND_KV_30DIFFSEED_NOIC", "HYPERVIS_SUBCYCLE", "RAD_BUG", "RAND-MT-DIFFSEED", "CPL_BUG", "EDG", "UW-SH", "GNU", "TSTEP_TYPE", "clm_init", "DUST", "FACTIC", "I15_03", "clm_albice00", "QSPLIT", "O0", "clm_veg", "RH-MIN-HIGH", "CONV-OCN", "CLDFRC_RHMINL", "clm_urban", "solar", "clm_default", "NU-P", "RH-MIN-LOW", "RAND_KV_DIFFSEED", "RAND-MT", "CONV_WATER"]
finished = [""]
for t in [x for x in testsUF if x not in finished]:
    try:
        testFiles = glob.glob(testDirUF + t + '/*.nc')
        test = sc.parallelize(testFiles, partitions)
        testTmp = test.map(gmsNetCDF).coalesce(partitions)
        testGms = testTmp.collect()
        ex, lbls = build_examples(testGms, gmsUF, 10000, aug=False)
        logistic = linear_model.RandomizedLogisticRegression(normalize=True, n_resampling=200, sample_fraction=0.75, selection_threshold=0.85, n_jobs=8)
        #logistic = linear_model.RandomizedLogisticRegression(n_resampling=200, selection_threshold=2.5)
        logistic.fit(ex, lbls)
        betas = logistic.scores_
        #betas = betas.reshape((len(testGms),len(varUF.value)))
        #avg_betas = np.mean(betas, axis=0)
        sort = np.argsort(-betas)#np.argsort(-avg_betas)#
        l = [varUF.value[i] for i in sort if betas[i] > 0.]#[varUF.value[i] for i in sort if avg_betas[i] > 0.]#
        print "Test set: ", t
        print "Number of important variables: ", len(l)
        print "Order of importance:"
        print l
        print "Alphabetical order:"
        print sorted(l)
        #bpd = distributions.main("/glade/p/work/milroy/sz750-ufect_T.nc", testGms, testDirUF + "/" + t + "/", "madeUpDir")
        bpd = distributions.main("/glade/p/work/milroy/sz750-ufect_T.nc", testDirUF + "/" + t + "/", testGms, "madeUpDir")
        print "Boxplot variables:"
        print bpd 
        rlr_vs_bp.main(l, bpd, [], [], t + "_noaug_resamp200_frac.75_norm_thresh.85_cast") 
        print "Experiment", t, "finished"
    except:
        pass

