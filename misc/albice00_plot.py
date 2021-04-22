import numpy as np
from collections import OrderedDict, defaultdict
import cPickle as pickle
import itertools
import seaborn as sns
import matplotlib.pyplot as plt

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

with open('albice00_sig.pickle', 'rb') as f:
  gms = pickle.load(f)

global_means = defaultdict(OrderedDict)

for var in gms['000']:
  global_means[0][var] = gms['000'][var]
  global_means[1][var] = gms['001'][var]

significance = OrderedDict()
for var in global_means[0]:
  sigfigs = np.zeros(global_means[0][var].shape[0], dtype=np.int)
  #precision = 16
  for i in xrange(0, 11):
    precision = 16
    while precision > 0:
      if RoundToSigFigs(global_means[0][var][i], precision) == RoundToSigFigs(global_means[1][var][i], precision):
        sigfigs[i] = precision
        break
      precision -= 1
  significance[var] = sigfigs

#with open('significance.pickle', 'rb') as f:
#  significance2 = pickle.load(f)

subcycled = ['QRL', 'QRS', 'AODDUST1', 'AODDUST3', 'AODVIS', 'BURDENBC', 'BURDENDUST', 'BURDENPOM', 
              'BURDENSEASALT', 'BURDENSO4', 'BURDENSOA', 'FLDS', 'FLNS', 'FLNSC', 'FLNT', 'FLNTC', 
              'FLUT', 'FLUTC', 'FSDS', 'FSDSC', 'FSNS', 'FSNSC', 'FSNT', 'FSNTC', 'FSNTOA', 'FSNTOAC']

for var in subcycled:
  significance[var][[1, 2, 4, 6, 8, 10]] = -1

mat = np.zeros((117, 11), dtype=np.int)
i = 0
for j in sorted(significance.keys()):
  mat[i, :] = significance[j]
  i += 1

masked_mat = np.ma.masked_where(mat == -1, mat)


fig, axn = plt.subplots(1, 2, figsize=(48., 42.))
cbar_ax = fig.add_axes([.91, .3, .03, .4])
colmap = plt.cm.get_cmap('YlOrRd_r', 15)
colmap.set_bad(color='black')
colmap.set_over(color="white")

data = [masked_mat.data[:58, :], masked_mat.data[58:, :]]
mask = [masked_mat.mask[:58, :], masked_mat.mask[58:, :]]
labels = [sorted(significance.keys())[:58], sorted(significance.keys())[58:]]

A = {}
for i, ax in enumerate(axn.flat):
    A[i] = sns.heatmap(data=data[i], xticklabels=range(0,11), 
                       yticklabels=labels[i], cmap=colmap, 
                       mask=mask[i], vmin=1, vmax=15, 
                       cbar_kws={"extend": "max", "drawedges": "true"}, ax=ax, cbar=i == 0, cbar_ax=None if i else cbar_ax)

axn[0].axhline(y=0, color='black',linewidth=2)
axn[0].axhline(y=masked_mat.data.shape[0]/2, color='black',linewidth=2)
axn[0].axvline(x=0, color='black',linewidth=2)
axn[0].axvline(x=masked_mat.data.shape[1], color='black',linewidth=2)

axn[1].axhline(y=0, color='black',linewidth=2)
axn[1].axhline(y=(masked_mat.data.shape[0]/2+1), color='black',linewidth=2)
axn[1].axvline(x=0, color='black',linewidth=2)
axn[1].axvline(x=masked_mat.data.shape[1], color='black',linewidth=2)

axn[0].set_xlabel('CESM time step', fontsize=58)
axn[0].xaxis.labelpad = 40
axn[1].set_xlabel('CESM time step', fontsize=58)
axn[1].xaxis.labelpad = 40
plt.sca(axn[0])
plt.yticks(rotation=0, fontsize=36)
plt.xticks(fontsize=52)
plt.sca(axn[1])
plt.yticks(rotation=0, fontsize=36)
plt.xticks(fontsize=52)

# removed = [2, 3, 4, 11, 12, 13, 14, 15, 16]
# i = 0
# for ytick in reversed(axn[0].get_yticklabels()):
#     if i in removed:
#         ytick.set_color('red')
#     i += 1

cbar = A[0].collections[0].colorbar
cbar.set_ticks([1, 3, 5, 7, 9, 11, 13, 15])
#cbar.set_ticks([1, 3, 5, 7, 9])
cbar.set_ticklabels(['1', '3', '5', '7', '9', '11', '13', '15'])#, pad=45)
#cbar.set_ticklabels(['1', '3', '5', '7', '9'])#, pad=45)
cbar.solids.set_edgecolor("black")
cax = plt.gcf().axes[-1]
cax.tick_params(labelsize=52)
tl = cax.get_yticklabels()
tl[0].set_verticalalignment('bottom')
tl[-1].set_verticalalignment('top')

fig.tight_layout(rect=[0, 0.05, .9, 1])
#fig.text(0.5, 0.025, 'CESM time step', ha='center', va='bottom', fontsize=58)

plt.savefig('albice00.pdf')

