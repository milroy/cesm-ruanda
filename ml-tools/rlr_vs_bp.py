import numpy as np
from collections import OrderedDict, defaultdict
import cPickle as pickle
import itertools
import seaborn as sns
import netCDF4 as nc
from matplotlib.colors import ListedColormap

def main(rlruf, bpuf, rlryr, bpyr, fname):
    baseuf = nc.Dataset('/glade/p/work/milroy/sz750-ufect_T.nc', 'r')
    varsuf = []

    baseyr= nc.Dataset('/glade/p/tdd/asap/verification/pca_ens_summary/cesm1_3_beta11/sz453.compilers-rand2_V6.nc', 'r')
    varsyr = []
    for i in xrange(baseuf.variables['vars'].shape[0]):
      varsuf.append(baseuf.variables['vars'][i,:].tostring().strip())

    varsuf = sorted(varsuf)
    for i in xrange(baseyr.variables['vars'].shape[0]):
      varsyr.append(baseyr.variables['vars'][i,:].tostring().strip())

    varsyr = sorted(varsyr)

    mat = np.zeros((len(varsyr), 2), dtype=np.int)
    if not (len(bpuf) == 0 and len(rlruf) == 0) and not (len(bpyr) == 0 and len(rlryr) == 0):
      i = 0
      for j in varsyr:
        if j in varsuf:
          if j in bpuf and j in rlruf:
            mat[i, 0] = 1
          elif j in bpuf:
            mat[i, 0] = 2
          elif j in rlruf:
            mat[i, 0] = 3
          else:
            mat[i, 0] = 4
        else:
          mat[i, 0] = 5

        if j in bpyr and j in rlryr:
          mat[i, 1] = 1
        elif j in bpyr:
          mat[i, 1] = 2
        elif j in rlryr:
          mat[i, 1] = 3
        else:
          mat[i, 1] = 4
        i += 1

    elif not (len(bpuf) == 0 or len(rlruf) == 0):
      i = 0
      for j in varsyr:
        if j in varsuf:
          if j in bpuf and j in rlruf:
            mat[i, 0] = 1
          elif j in bpuf:
            mat[i, 0] = 2
          elif j in rlruf:
            mat[i, 0] = 3
          else:
            mat[i, 0] = 4
        else:
          mat[i, 0] = 5
        mat[i, 1] = 5
        i += 1

    elif not (len(bpyr) == 0 or len(rlryr) == 0):
      i = 0
      for j in varsyr:
        if j in bpyr and j in rlryr:
          mat[i, 1] = 1
        elif j in bpyr:
          mat[i, 1] = 2
        elif j in rlryr:
          mat[i, 1] = 3
        else:
          mat[i, 1] = 4
        mar[i, 0] = 5
        i += 1

    fig, axn = sns.plt.subplots(1, 2, figsize=(48., 42.))
    cbar_ax = fig.add_axes([.91, .3, .03, .4])
    colmap = ListedColormap(['purple', 'blue', 'red', 'white', 'black'])

    data = [mat[:len(varsyr)/2, :], mat[len(varsyr)/2:, :]]
    mask = [mat[:len(varsyr)/2, :], mat[len(varsyr)/2:, :]]
    labels = [varsyr[:len(varsyr)/2], varsyr[len(varsyr)/2:]]

    A = {}
    for i, ax in enumerate(axn.flat):
        A[i] = sns.heatmap(data=data[i], xticklabels=['UF', 'YR'], 
                           yticklabels=labels[i], cmap=colmap, 
                           vmin=1, vmax=5, 
                           ax=ax, cbar=i == 0, cbar_ax=None if i else cbar_ax)

    axn[0].axhline(y=0, color='black',linewidth=2)
    axn[0].axhline(y=mat.shape[0]/2, color='black',linewidth=2)
    axn[0].axvline(x=0, color='black',linewidth=2)
    axn[0].axvline(x=mat.shape[1], color='black',linewidth=2)

    axn[1].axhline(y=0, color='black',linewidth=2)
    axn[1].axhline(y=(mat.shape[0]/2+1), color='black',linewidth=2)
    axn[1].axvline(x=0, color='black',linewidth=2)
    axn[1].axvline(x=mat.shape[1], color='black',linewidth=2)

    axn[0].set_xlabel('Vars 1-58', fontsize=58)
    axn[0].xaxis.labelpad = 40
    axn[1].set_xlabel('Vars 59-117', fontsize=58)
    axn[1].xaxis.labelpad = 40
    sns.plt.sca(axn[0])
    sns.plt.yticks(rotation=0, fontsize=36)
    sns.plt.xticks(fontsize=52)
    sns.plt.sca(axn[1])
    sns.plt.yticks(rotation=0, fontsize=36)
    sns.plt.xticks(fontsize=52)

    cbar = A[0].collections[0].colorbar
    cbar.set_ticks([1, 2, 3, 4, 5])
    cbar.set_ticklabels(['Both', 'BP', 'RLR', 'Neither', 'Not computed'])#, pad=45)
    #cbar.solids.set_edgecolor("black")
    cax = sns.plt.gcf().axes[-1]
    cax.tick_params(labelsize=52)
    tl = cax.get_yticklabels()
    tl[0].set_verticalalignment('bottom')
    tl[-1].set_verticalalignment('top')

    fig.tight_layout(rect=[0, 0.05, .9, 1])
    #fig.text(0.5, 0.025, 'CESM time step', ha='center', va='bottom', fontsize=58)

    sns.plt.savefig(fname + ".pdf")
    sns.plt.clf()


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Dist script')
  parser.add_argument('rlruf', type=list, default='')
  parser.add_argument('bpuf', type=list, default='')
  parser.add_argument('rlryr', type=list, default='')
  parser.add_argument('bpyr', type=list, default='')
  parser.add_argument('fname', type=str, default='')

  args = parser.parse_args()

  main(args.rlruf, args.bpuf, args.rlryr, args.bpyr, args.fname)


