import numpy as np
import glob
import os
import cPickle as pickle
import itertools
import Nio as nio
from pyspark.sql import SparkSession
spark = SparkSession.builder.getOrCreate()
sc = spark.sparkContext
import pyspark
from pyspark import SparkConf, SparkContext
import numpy as np
from collections import OrderedDict, defaultdict
import re

def gmsNetCDF(part, type):
    f = nio.open_file(part, 'r')
    filenumber = re.findall(r"\D(\d{3})\D", part)[0]
    area = f.variables["area"]
    area_wgt = area[:].astype(np.float64)
    total = np.sum(area_wgt)
    area_wgt /= total
    nlev = f.dimensions["lev"]
    k = 0
    if type == "uf":
        global_means_run = np.zeros((len(varUF.value), f.variables[varUF.value[0]].shape[0]), dtype=np.float64)
        itervar = varUF.value
    elif type == "yr":
        global_means_run = np.zeros((len(varYr.value), f.variables[varYr.value[0]].shape[0]), dtype=np.float64)
        itervar = varYr.value
    for vari in itervar:
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
    return filenumber, global_means_run


partitions = 1
baseYr = nio.open_file('/glade/p/tdd/asap/verification/pca_ens_summary/cesm1_3_beta11/sz453.compilers-rand2_V6.nc', 'r')
varssYr = []
for i in xrange(baseYr.variables['vars'].shape[0]):
    varssYr.append(baseYr.variables['vars'][i,:].tostring().strip())

varYr = sc.broadcast(varssYr)

testFiles = glob.glob('/glade/scratch/milroy/ufect.hydro_baseflow' + '/*.nc')
test = sc.parallelize(testFiles, partitions)
testTmp = test.map(lambda file: gmsNetCDF(file, "yr")).coalesce(partitions)
testGms = testTmp.collect()

test = {}
for i in testGms:
  sim = OrderedDict()
  for j in xrange(i[1].shape[0]):
    sim[varYr.value[j]] = i[1][j]
  test[i[0]] = sim

__logBase10of2 = 3.010299956639811952137388947244930267681898814621085413104274611e-1
def RoundToSigFigs(x, sigfigs):
    if not ( type(sigfigs) is int or np.issubdtype(sigfigs, np.integer)):
        raise TypeError( "RoundToSigFigs: sigfigs must be an integer." )
    if not np.all(np.isreal( x )):
        raise TypeError( "RoundToSigFigs: all x must be real." )
    if sigfigs <= 0:
        raise ValueError( "RoundtoSigFigs: sigfigs must be positive." )
    mantissas, binaryExponents = np.frexp(x)
    decimalExponents = __logBase10of2 * binaryExponents
    intParts = np.floor(decimalExponents)
    mantissas *= np.power(10.0, (decimalExponents - intParts))
    return np.around(mantissas, decimals=sigfigs - 1) * np.power(10.0, intParts)

global_means_test = test
tstep = 0
l = {}
for pert in global_means:
  significance = OrderedDict()
  for var in global_means[pert]:
    if not global_means[pert][var][tstep] == global_means_test[pert][var][tstep]:
      precision = 16
      while precision > 0:
        if RoundToSigFigs(global_means[pert][var][tstep], precision) == RoundToSigFigs(global_means_test[pert][var][tstep], precision):
          break
        precision -= 1
      significance[var] = precision
  l[pert] = significance

