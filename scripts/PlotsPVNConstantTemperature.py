import matplotlib.pyplot as plt
import numpy as np
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-d', '--result-directory', dest='Result_Directory', type=str)
parser.add_argument('-o', '--output-filename', dest='OutFileName', type=str)
args = parser.parse_args()

plt.style.use('/Users/Nappo/work/local/mplstyles/root_style.mplstyle')

for f in sorted(os.listdir(args.Result_Directory)):
  if(f.endswith('.txt') and f.startswith('results')):
    temp = float(f.replace('results','').replace('.txt',''))
    data = np.loadtxt(args.Result_Directory + '/' + f ,
                      delimiter=',').transpose()
    print(data[2])
    plt.errorbar(data[0],data[1],yerr=data[2],
      marker='o',label=r'$T={0:.1f}$'.format(temp), markersize=3)
plt.ylabel(r'$P$')
plt.xlabel(r'$V/N$')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig(args.OutFileName + '.pdf')