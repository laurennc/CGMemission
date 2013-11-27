from lauren import *

#what types of plots do I want to make.....

#recreate Bertone 2010 Figure 1 where I'm plotting cloudy results for different densities

def plot_Cloudy_output(inputfile,outputfile):
	cldata = np.genfromtxt(inputfile,skip_header=??,dtype=[('Temperature',float),('Lya',float),('Ha',float),('C4a',float),('C4b',float),('O6a',float),('O6b',float)])
	cldata = cldata.view(np.recarray)
	
	data = np.zeros((len(cldata['Temperature']),5))
	data[:,0] = cldata['Temperature']
	data[:,1] = cldata['Lya']
	data[:,2] = cldata['Ha']
	data[:,3] = np.log10(10.**(cldata['C4a'])+10.**(data['C4b']))
	data[:,4] = np.log10(10.**(cldata['O6a'])+10.**(data['06b']))

	labels = ['Temperature','Lya','Ha','CIV','OVI']
	i = 1
	while (i < 5):
		plt.plot(data[:,0],data[:,i],label=labels[i],linewidth=3)
		i = i + 1
	plt.xlim()
	plt.ylim()
	plt.legend()
	plt.xlabel(r"")
	plt.ylabel(r"")
	plt.savefig(outputfile)

#make comparison for the cloudy runs I DO have of the scaled values and the actual Cloudy values to see how they compare
#NOTE: I'm scaling simply by the metallicity and not by a assumed abundance based off of that metallicity



