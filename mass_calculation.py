import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from yt.mods import *
from lauren import distance_from_center

fn = "/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
pf = load(fn, file_style="%s.grid.cpu%%04i")


#val,pos = pf.h.find_max('Density')
pos = [ 0.40328598,  0.47176743,  0.46131516]

pos2 = [ 0.40320258,  0.4718025 ,  0.46123817]
rad = 316.8648 ##kpc

data = pf.h.sphere(pos2,rad/pf['kpc'])
dists = distance_from_center(data['x'],data['y'],data['z'],pos)*pf['kpc']

plt.plot(dists,np.log10(data['Density']),'b.',alpha=0.15)
plt.xlabel('Radius [kpc]')
plt.xlim(0,30)
plt.ylabel('log(Density)')
plt.savefig('dens_radius_fordisk_denscent.png')
plt.close()

dothis = False
if dothis:
	idx = np.where(np.log10(data['Temperature']) > 6.0)[0]
	hotmass = np.log10(np.sum(data['CellMassMsun'][idx]))

	idx = np.where((np.log10(data['Temperature']) => 5.0) & (np.log10(data['Temperature']) <= 6.0 ) )[0]
	warmmass = np.log10(np.sum(data['CellMassMsun'][idx]))

	idx = np.where(np.log10(data['Temperature']) < 5.0 )[0]
	coldmass = np.log10(np.sum(data['CellMassMsun'][idx]))

	



