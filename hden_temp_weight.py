import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import cPickle
from yt.mods import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MultipleLocator


emis = cPickle.load(open('cloudywerk.cpkl','rb'))

weight_col = 'SiIV'
weight_basic = 'SiIV_Density'

fileout = 'histdens_z0_g1q1_'+weight_basic+'.png'

hden_file = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/basic/z0/g1q1/frbx_1kpc_500kpc_z0_H_NumberDensity_'+weight_basic+'.cpkl'
temp_file = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/basic/z0/g1q1/frbx_1kpc_500kpc_z0_Temperature_'+weight_basic+'.cpkl'
column_file = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g1q1/frbx_1kpc_500kpc_z0_'+weight_col+'dens.cpkl'

hden = cPickle.load(open(hden_file,'rb'))
temp = cPickle.load(open(temp_file,'rb'))
coldens = cPickle.load(open(column_file,'rb'))

hden,temp,coldens = np.log10(hden.flat),np.log10(temp.flat),np.log10(coldens.flat)

hist, xedges, yedges  = np.histogram2d(hden,temp,weights=coldens,range=[[-6,1],[4, 6]],bins=100)

xbins = np.digitize(hden,xedges[1:-1])
ybins = np.digitize(temp,yedges[1:-1])

totvals = np.zeros((len(xedges[1:-1]),len(yedges[1:-1])))
numvals = np.zeros((len(xedges[1:-1]),len(yedges[1:-1])))
maxvals = np.zeros((len(xedges[1:-1]),len(yedges[1:-1])))

for i in range(len(coldens)):
    totvals[xbins[i],ybins[i]] += 10.**coldens[i]
    numvals[xbins[i],ybins[i]] += 1
    if 10.**coldens[i] > maxvals[xbins[i],ybins[i]]:
        maxvals[xbins[i],ybins[i]] = 10.**coldens[i]

avgvals = totvals/numvals

idx = np.where(totvals == 0)
avgvals[idx] = 0

fig,ax = plt.subplots(1,3)
fig.set_size_inches(12,4)
ax = ax.flat

im1 = ax[0].imshow(np.log10(hist.T),extent=[-6,1,4,6],interpolation='nearest',origin='lower',cmap='Blues')
im2 = ax[1].imshow(np.log10(avgvals.T),extent=[-6,1,4,6],interpolation='nearest',origin='lower',vmin=9.,vmax=16.)
im3 = ax[2].imshow(np.log10(maxvals.T),extent=[-6,1,4,6],interpolation='nearest',origin='lower',vmin=9.,vmax=16.)

x0,x1 = ax[2].get_xlim()
y0,y1 = ax[2].get_ylim()
titles = ['Number','Avg Column Density','Max Column Density']

for j in range(len(ax)):
   ax[j].set_aspect((x1-x0)/(y1-y0))
   ax[j].set_xlim(-6,1)
   ax[j].set_ylim(3.5,6)
   ax[j].set_title(titles[j])
  

ax[0].plot(emis['hden'],emis['T'],'s',color='HotPink',alpha=0.75,ms=3.5)


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im3, cax=cbar_ax)#,shrink=0.6,ticks=MultipleLocator(0.5), format="%.2f")


#plt.set_cmap('spectral')
plt.savefig(fileout)
plt.close()

