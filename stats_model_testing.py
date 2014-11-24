import matplotlib
matplotlib.use('Agg')
import cPickle
import numpy as np
import matplotlib.pyplot as plt


def data_number_vs_rperp(r_data,r_tol):
	xL = np.linspace(-160,160,320)
	xL,yL = np.meshgrid(xL,xL)
	radii = abs(xL+1j*yL)
	num_pts = []
	
        for r in r_data:
                idr = np.where((radii >= r-r_tol) & (radii <= r+r_tol))
                rad = radii[idr]
                num_pts = np.append(num_pts,len(rad))

	return num_pts

galaxy_props = cPickle.load(open('werk_galaxy_properties.cpkl','rb'))

num_pts = data_number_vs_rperp(galaxy_props['Rperp'],2.5)

plt.plot(galaxy_props['Rperp'],num_pts,'bo')
plt.xlabel('Rperp (kpc)')
plt.ylabel('Number of Model Points')
plt.savefig('available_model_points_2kpc.png')
plt.close()


