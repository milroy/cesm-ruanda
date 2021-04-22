#!/usr/bin/python

import numpy as np
from collections import OrderedDict, defaultdict
import Nio as nio
import cPickle as pickle
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
import gc

base = nio.open_file('/home/milroy/NCAR/sz750-ufect_T.nc', 'r')

vars_s = []

for i in xrange(base.variables['vars'].shape[0]):
    vars_s.append(base.variables['vars'][i,:].tostring().strip())

with open('/home/milroy/NCAR/var_tests/ens_10ts.pickle', 'rb') as f:
    global_means = pickle.load(f)

with open('/home/milroy/NCAR/var_tests/albice00_10ts.pickle', 'rb') as f:
    global_means_test = pickle.load(f)

boolean = []
for i in xrange(11):
  if i%2 == 0 and not i == 0 or i == 1:
    boolean.append(True)
  else:
    boolean.append(False)

subcycled = ['FLDS', 'FLNS', 'FLNSC', 'FLNT', 'FLNTC', 'FLUT', 'FLUTC',
              'FSDS', 'FSDSC', 'FSNS', 'FSNSC', 'FSNT', 'FSNTC', 'FSNTOA',
              'FSNTOAC', 'QRL', 'QRS']

diff_scaled = OrderedDict()
for sim in sorted(global_means.keys()):
    print sim
    diffvar = OrderedDict()
    for var in global_means[sim]:
        if var in subcycled:
            #print var
            ma = (global_means_test[sim][var] - global_means[sim][var])/(global_means[sim][var])
            diffvar[var] = np.ma.array(ma, mask=boolean)
        else:
            diffvar[var] = (global_means_test[sim][var] - global_means[sim][var])/(global_means[sim][var])
    diff_scaled[sim] = diffvar


# diff_scaled = OrderedDict()
# for sim in sorted(global_means.keys()):
#     print sim
#     diffvar = OrderedDict()
#     for var in global_means[sim]:
#         diffvar[var] = (global_means_test[sim][var] - global_means[sim][var])/(global_means[sim][var])
#     diff_scaled[sim] = diffvar

k = 1
for vary in vars_s:
    palette = itertools.cycle(sns.color_palette(sns.cubehelix_palette(30, start=.5, rot=-.75)))
    plt.figure(figsize=(32.,18.))
    sns.set(font_scale=2)
    sns.set_style("darkgrid")
    plt.gcf().suptitle(vary, fontsize=45)
    if vary in subcycled:
        for sim in diff_scaled:
            plot = sns.tsplot(data=diff_scaled[sim][vary][~diff_scaled[sim][vary].mask], color=next(palette), linewidth=3)
    else:
        for sim in diff_scaled:
            plot = sns.tsplot(data=diff_scaled[sim][vary], color=next(palette), linewidth=3)
    name = vary + '.pdf'
    plt.tight_layout()
    plt.savefig('var-ts/' + name)
    plt.savefig('var-ts/' + str(k) + '.png')
    plt.clf()
    plt.close()
    gc.collect()
    k += 1
