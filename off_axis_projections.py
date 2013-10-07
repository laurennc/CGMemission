import matplotlib as mpl
mpl.use('agg')
from yt.mods import *
import numpy as np
import triangle
from matplotlib.artist import *

fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

#This is the position of maximum density
pos = [0.40328598,0.47176743,0.46131516]
#This radius is chosen to visually match the plots in Ximena's paper
rad = 108.0/pf['kpc']
data = pf.h.sphere(pos,rad)

#bah want the gas density....
densities = np.log10(data['Density'])
temperatures = np.log10(data['Temperature'])

idx = np.where((temperatures > 3.8) & (temperatures < 4.25) )[0]
densities = densities[idx]
temperatures = temperatures[idx]
xpos = data['x'][idx]
ypos = data['y'][idx]
zpos = data['z'][idx]

idx2 = np.where(densities > -24)[0]
densities = densities[idx2]
temperatures = temperatures[idx2]
xpos = xpos[idx2]
ypos = ypos[idx2]
zpos = zpos[idx2]

datain = np.zeros((len(densities),2))
datain[:,0] = densities
datain[:,1] = temperatures

#fileout = 'dens_temp3.png'
#triangle.corner(datain,labels=['log(Density)','log(Temperature)']).savefig(fileout)

fileout = 'disk_positions2.png'
datain2 = np.zeros((len(densities),3))
datain2[:,0] = xpos#*pf['kpc']
datain2[:,1] = ypos#*pf['kpc']
datain2[:,2] = zpos#*pf['kpc']
triangle.corner(datain2,labels=['x','y','z']).savefig(fileout)


