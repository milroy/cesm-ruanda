#!/usr/bin/python

import numpy as np
from collections import OrderedDict
from collections import defaultdict
import cPickle as pickle
import itertools
import seaborn as sns
import pandas as pd
import Nio as nio
import glob
import re
import os
import argparse
import pdb

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
    #print filename
    filenumber = re.findall(r"\D(\d{3})\D", filename)[0]
    fn = int(filenumber)
    #print fn
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
      #print sim
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
      #print sim
      diffvar = OrderedDict()
      for var in ens_global_means[sim]:
        diffvar[var] = (ens_global_means[sim][var][tstep] - ens_global_means[0][var][tstep])/(ens_global_means[0][var][tstep])
      diff_scaled_ens[sim] = diffvar

  ens = defaultdict(OrderedDict)
  for sim in diff_scaled_ens:
    for var in diff_scaled_ens[sim]:
      ens[var][sim] = diff_scaled_ens[sim][var]

  return (ens, exp)


def plot(ens, exp, enslabel, explabel, path):
  colors = sns.color_palette("colorblind")

#  for case in ['case1', 'case2', 'case3']:
#    if not os.path.exists(path + '/' + case):
#      os.makedirs(path + '/' + case)
  case1 = []
#  case2 = []
#  case3 = []
  for vary in ens:
    #print exp[vary]
    exp_l = [exp[vary][x] for x in exp[vary]]
    ens_l = [ens[vary][x] for x in ens[vary]]
    exp_median = np.median(exp_l)
    ens_median = np.median(ens_l)
    iqr_exp = np.percentile(exp_l, [75, 25])
    iqr_ens = np.percentile(ens_l, [75, 25])

#    data = pd.DataFrame.from_dict(OrderedDict([(enslabel, ens_l), (explabel, exp_l)]))
#    sns.set_style("darkgrid")
#    sns.boxplot(data=data, palette=colors, whis=[2.5, 97.5])
#    sns.plt.suptitle(vary)

    if iqr_exp[1] > iqr_ens[0] or iqr_ens[1] > iqr_exp[0]:
      case1.append(vary)
    elif (exp_median > iqr_ens[0] or exp_median < iqr_ens[1]) and (ens_median > iqr_exp[0] or ens_median < iqr_exp[1]):
#      case2.append(vary)
      pass
    else:
#      case3.append(vary)
      pass
  return case1

def main(sumFile, testDir, testGms, outDir, testFile, tstep):
  sumFile = str(sumFile)
  testDir = str(testDir)
  outDir = str(outDir)
  testFile = str(testFile)

  if testDir:
    numtestFiles = len(os.listdir(testDir))

  elif testFile:


  base = nio.open_file(sumFile, 'r')

  path_exp = testDir
  exp_gms, ens_gms = get_gms(base, testGms, path_exp)
  scaled_ens, scaled_exp = scale(ens_gms, exp_gms, numtestFiles, tstep)

#  pdb.set_trace()
  return plot(scaled_ens, scaled_exp, 'ensemble', 'experiment', outDir)


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Dist script')
  parser.add_argument('sumFile', type=str, default='')
  parser.add_argument('testDir', type=str, default='')
  parser.add_argument('testFile', type=str, default='')
  parser.add_argument('testGms', type=list)
  parser.add_argument('outDir', type=str, default='')
  parser.add_argument('tstep', type=int, default=0)


  args = parser.parse_args()

  main(args.sumFile, args.testDir, args.testGms, args.outDir, args.testFile, args.tstep)


