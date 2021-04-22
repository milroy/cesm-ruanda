#!/usr/bin/python

import numpy as np
from collections import OrderedDict
from collections import defaultdict
import cPickle as pickle
import itertools
import seaborn as sns
import pandas as pd
import numpy as np

with open('albice00_60_tstep.pickle', 'rb') as f:
  global_means = pickle.load(f)

with open('60_tstep.pickle', 'rb') as f:
  control = pickle.load(f)

boolean = []
for i in xrange(61):
  
  if i%2 == 0 and not i == 0 or i == 1:
    boolean.append(True)
  else:
    boolean.append(False)

subcycled = ['FLDS', 'FLNS', 'FLNSC', 'FLNT', 'FLNTC', 'FLUT', 'FLUTC',
              'FSDS', 'FSDSC', 'FSNS', 'FSNSC', 'FSNT', 'FSNTC', 'FSNTOA', 
              'FSNTOAC', 'QRL', 'QRS']

diff_scaled = OrderedDict()
for sim in global_means:
    print sim
    diffvar = OrderedDict()
    for var in global_means[sim]:
        if var in subcycled:
            #print var
            ma = (global_means[sim][var] - control['000'][var])/(control['000'][var])
            diffvar[var] = np.ma.array(ma, mask=boolean)
        else:
            diffvar[var] = (global_means[sim][var] - control['000'][var])/(control['000'][var])
    diff_scaled[sim] = diffvar

tmp1 = defaultdict(OrderedDict)
for sim in diff_scaled:
  for var in diff_scaled[sim]:
      tmp1[var][sim] = diff_scaled[sim][var]

diff_scaledc = OrderedDict()
for sim in control:
    diffvar = OrderedDict()
    for var in control[sim]:
        if var in subcycled:
            ma = (control[sim][var][1] - control['000'][var][1])/(control['000'][var][1])
            diffvar[var] = np.ma.array(ma, mask=boolean)
        else:
            diffvar[var] = (control[sim][var][1] - control['000'][var][1])/(control['000'][var][1])
    diff_scaledc[sim] = diffvar

tmp2 = defaultdict(OrderedDict)
for sim in diff_scaledc:
  for var in diff_scaledc[sim]:
      tmp2[var][sim] = diff_scaledc[sim][var]

for vary in global_means['000']:
    sns.plt.figure(figsize=(32.,24.))
    colors = sns.color_palette("colorblind")
    sns.set(font_scale=2)
    sns.set_style("darkgrid")
    palette = itertools.cycle(colors)
    if vary in subcycled:
        for sim in diff_scaled:
            if sim != '000':
                print sim
                c = next(palette)
                sns.tsplot(data=diff_scaled[sim][vary][~diff_scaled[sim][vary].mask], color="#e74c3c")
                sns.tsplot(data=diff_scaledc[sim][vary][~diff_scaledc[sim][vary].mask], color="#34495e")
    else:
        for sim in diff_scaled:
            if sim != '000':
                print sim
                c = next(palette)
                sns.tsplot(data=diff_scaled[sim][vary], color="#e74c3c")
                sns.tsplot(data=diff_scaledc[sim][vary], color="#34495e")        
    sns.plt.suptitle(vary)
    sns.plt.savefig('albice_60ts/' + vary + '.pdf')
    sns.plt.clf()
