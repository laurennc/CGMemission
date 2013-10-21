import matplotlib as mpl
mpl.use('agg')
from lauren import *
import matplotlib.pyplot as plt
import triangle as triangle

fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"            
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data
pos = [0.40328598,0.47176743,0.46131516]

#David would like me to make different cuts on Temperature and Volume to see what effects we're seeing --> seemed to like the T,HI,EM,Ha plots so continue those

#what volume cuts are interesting? 20,50,80,100
#what temperature cuts are interesting? 4.5, 5.0, 5.5, 6.0
#cell be cell is much much easier than projections...

radii = [20.,50.,80.,100.]
temps = [4.5,5.0,5.5,6.0]
fields = ['Temp','HI Col Dens','EM', 'Ha Arc']
fieldlims = [(4,8),(8,18),(-24,-14),(-46,-36)]

radcount = 0

for rad in radii:
	data = pf.h.sphere(pos,rad/pf['kpc'])	
	count = 0
	fig, axes = plt.subplots(4,6)
	fig.set_size_inches(10,10)
	
	#want to have only annuli
	if rad > radii[0]:
		dists = distance_from_center(data['x'],data['y'],data['z'],pos)
		idD = np.where((dists > radii[radcount-1]/pf['kpc']) & (dists <= radii[radcount]/pf['kpc'] ))[0]

	for temp in temps:
		idx = np.where(np.log10(data['Temperature'][idD]) > temp)[0]
		if rad > radii[0]:
			idx = np.intersect1d(idD,idx)
		datain = np.zeros((len(idx),4))
		datain[:,1] = np.log10(data['HI_Number_Density'][idx]*data['dx'][idx]*pf['cm'])
		datain[:,2] = np.log10(data['EmissionMeasureCM'][idx]*data['dx'][idx]*pf['pc'])
		datain[:,3] = np.log10(data['HAlphaEmissionArc'][idx]*data['dx'][idx]*pf['cm'])
		datain[:,0] = np.log10(data['Temperature'][idx])

		#ok so this is fine but how do I make the numbers fit correctly?
		otro = 0
		j = 0
		while j < 3:
			i = j + 1
			while i < 4:
				print 'Plot axes: ',count,otro
				print 'Data axes: ',j, i
				ax = axes[count,otro]
				triangle.hist2d(datain[:,j],datain[:,i],ax=ax) 
				ax.locator_params(axis='x',nbins=4)
				ax.set_xlim(fieldlims[j])
				ax.set_ylim(fieldlims[i])
				if count == 0:
					ax.set_title(fields[i])
				if count== 3:
					ax.set_xlabel(fields[j])
				if otro == 0:
					ax.set_ylabel('Temp < '+str(temp))
				#ax.locator_params(axis='x',tight=True,nbins=5)
				otro = otro + 1
				i = i + 1
			j = j + 1
		count = count + 1

	#plt.suptitle('Data Sphere with Radius of '+str(rad)+' kpc')
	radcount = radcount + 1
	plt.tight_layout(w_pad=1.2,h_pad=1.2)
	fileout = ('varyTemp_rad='+str(rad)+'kpc.png')
	plt.savefig(fileout)
	plt.close()

#python can merge plots into a pdf?

