import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from plotting_routines import *



ions = ['H','C','O']
ionnum = [1,3,6]
fig,ax = plt.subplots(1,3)
ax = ax.flatten()

for i in range(len(ions)):

        pattbeg = '/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/Ions/grid_galquas'
        inputfiles = ['/g01q01/bert_ion_run3_','/g1q1/bert_ion_run7_','/g10q10/bert_ion_run9_']
        txtDens = ['n=-5','n=-3','n=-2']
        colors = ['r','g','b']

        for j in range(len(inputfiles)):
                inputfile = pattbeg + inputfiles[j] + ions[i] + '.dat'
		print inputfile
                plot_specific_ion_fraction(inputfile,ax[i],ionnum[i],colors[j],labelOn=False,txt='',plotDashed=False,plotDotted=False)

plt.savefig('constant_gamma_n.png')
plt.close()


