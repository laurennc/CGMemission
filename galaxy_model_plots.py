import matplotlib
matplotlib.use('Agg')
from Galaxy import *
import cPickle
import numpy as np
import matplotlib.pyplot as plt

galaxies = cPickle.load(open('galaxy_models.cpkl','rb'))

galaxy_props = cPickle.load(open('werk_galaxy_properties.cpkl','rb'))

q = []
mass = []
rperp = []
for galaxy in galaxies:
	ID = galaxy.ID
	idx = np.where(galaxy_props['ID'] == ID)[0]
	if len(idx)!= 0:
		q = np.append(q,galaxy.best_q)
		mass = np.append(mass,galaxy_props['mStar'][idx])
		rperp = np.append(rperp,galaxy_props['Rperp'][idx])


fig, ax = plt.subplots(1,2,sharey=True)

print mass.shape,q.shape,rperp.shape
ax[0].plot(mass,np.log10(q),'bo')
ax[0].set_xlabel('Mstar')
ax[0].set_ylabel('q')

ax[1].plot(rperp,np.log10(q),'gs')
ax[1].set_xlabel('Rperp (kpc)')

plt.savefig('galprops_vs_q.png')
plt.close()




