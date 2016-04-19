import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import cPickle
import itertools
from werk_galaxies import galaxies
from yt.mods import *


def plot_radial_profiles(sphere,color,z_label):
	rad_profile = BinnedProfile1D(sphere, 100, "Radiuskpc", 0.0, 500., log_space=False)
	rad_profile.add_fields(["Density","Temperature","Metallicity"])
	i = 0
	while i < len(fields):
		#j,ixs = 0,[]
		#while j < len(rad_profile["Radiuskpc"])-1:
			#if np.log10(rad_profile[fields[i]][j+1] ) < (np.log10(rad_profile[fields[i]][j]) +0.05):
			#	ixs = np.append(ixs,int(j+1))
			
			#j = j + 1
	
		#print i,len(ixs)
		#print ixs[0:10]	
		#ixs = [int(x) for x in ixs]
		#print ixs[0:10]
		if i == (len(fields)-1):
			print 'Here I am'
			ax[i].plot(rad_profile["Radiuskpc"],np.log10(rad_profile[fields[i]]),linewidth=1.5,color=color,label=z_label)
		else:
			ax[i].plot(rad_profile["Radiuskpc"],np.log10(rad_profile[fields[i]]),linewidth=1.5,color=color)
		i = i + 1

	return
   


fields = ['Density','Temperature','Metallicity']
units = ['[g cm^-3]','[K]','[Z/Z_sun]']
fig,ax = plt.subplots(1,3)
fig.set_size_inches(12,4)
## 0, 0.2, 0.5, 1.0
colors = ['#e7298a','#d95f02','#1b9e77','#7570b3','#e7298a']

fn_beg = "/u/10/l/lnc2115/vega/data/Ryan/"
fns = ["r0058_l10/redshift0058","r0054/redshift0054",'r0048/redshift0048',"r0038/redshift0038"]
poss = [[0.40320258,0.4718025,0.46123817],[0.39898226,0.46904468,0.46830814],[0.3932969,0.46559631,0.47849673],[0.38580719,0.46065774,0.49183819]]
rads = [320.,320.,320.,320.]
zs = ['0','0.2','0.5','1']

for i in range(len(rads)):
	fn = fn_beg + fns[i]
	pf = load(fn, file_style="%s.grid.cpu%%04i")
	sphere = pf.h.sphere(poss[i],rads[i]/pf['kpc'])

	plot_radial_profiles(sphere,colors[i],zs[i])


for i in range(len(fields)):
	_ = ax[i].set_xlabel('Radius [kpc]')
	_ = ax[i].set_ylabel(fields[i]+' '+units[i])

ax[len(fields)-1].legend()	

plt.tight_layout()
plt.savefig('radial_profiles_referee_original.pdf')
plt.close()


