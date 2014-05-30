from lauren import *

#what types of plots do I want to make.....

#recreate Bertone 2010 Figure 1 where I'm plotting cloudy results for different densities

def readin_Cloudy_data(inputfile):
	cldata = np.genfromtxt(inputfile,skip_header=13,dtype=[('Temperature',float),('Lya',float),('Ha',float),('C4a',float),('C4b',float),('O6a',float),('O6b',float),('C3a',float),('C3b',float),('C3c',float)])
        cldata = cldata.view(np.recarray)

        data = np.zeros((len(cldata['Temperature']),7))
        data[:,0] = cldata['Temperature']
        data[:,1] = cldata['Lya']
        data[:,2] = cldata['Ha']
        data[:,3] = np.log10(10.**(cldata['C4a'])+10.**(cldata['C4b']))
        data[:,4] = np.log10(10.**(cldata['O6a'])+10.**(cldata['O6b']))
	data[:,5] = cldata['C3a']
	data[:,6] = np.log10(10.**(cldata['C3b'])+10.**(cldata['C3c']))
	return data

def plot_Cloudy_output(inputfile,axes,scale=False,scale_value=-0.3):#outputfile):
	data = readin_Cloudy_data(inputfile)

	labels = ['Temperature','Lya','Ha','CIV','OVI','CIII 977','CIII']
	i = 1
	while (i < 7):
		if scale==True:
			#print 'in scale loop'
			#print scale_value
			data[:,i] = np.log10(scale_by_metallicity(10.**(data[:,i]),-0.3,scale_value))
		axes.plot(data[:,0],data[:,i],label=labels[i],linewidth=3)
		i = i + 1
	axes.set_autoscale_on(False)
	axes.set_aspect('equal')
	axes.set_xlim(3,8)
	axes.set_ylim(-28,-20)#-22)
	#axes.legend(loc="lower right")
	axes.set_xlabel(r"log(Temperature)")
	axes.set_ylabel("\epsilon / n_{H}^2")
	#plt.savefig(outputfile)
	return

def plot_Cloudy_loop(inputfileS,outputfile,xlen,ylen):
	fig,ax = plt.subplots(ylen,xlen,sharex=True,sharey=True)
	fig.set_size_inches(11,4.5)
	fig.subplots_adjust(hspace=0.1,wspace=0.1)
	i = 0
	count = 0
	while i < xlen:
		#j = 0
		#while j < ylen:
			#plot_Cloudy_output(inputfileS[count],ax[i,j])
			#j = j + 1
			#count = count + 1
		plot_Cloudy_output(inputfileS[i],ax[i],scale_value=0.0)
		i = i + 1
	ax[0].text(7,-22.5,'n = 1')
	ax[1].text(7,-22.5,'n = -3')
	ax[2].text(7,-22.5,'n = -6')
	#ax[xlen-1].legend(loc="lower right")
	ax[0].legend(loc="lower right")
	plt.savefig(outputfile)
	return 

#make comparison for the cloudy runs I DO have of the scaled values and the actual Cloudy values to see how they compare
#NOTE: I'm scaling simply by the metallicity and not by a assumed abundance based off of that metallicity
def plot_Cloudy_scaling_compare(originalfile1,newmetalfile2,scaleTo,outputfile):
	fig,ax = plt.subplots(1,3)
	fig.set_size_inches(10,5)
	fig.subplots_adjust(wspace=0.3)
	plot_Cloudy_output(newmetalfile2,ax[0])
	plot_Cloudy_output(originalfile1,ax[1],scale=True,scale_value=0)
	
	data1 = readin_Cloudy_data(originalfile1)
	data2 = readin_Cloudy_data(newmetalfile2)
	
	labels = ['Temperature','Lya','Ha','CIV','OVI']
        i = 1
        while (i < 5):
		data1[:,i] = np.log10(scale_by_metallicity(10.**data1[:,i],-0.3,0))
                ax[2].plot(data1[:,0],10.**(data1[:,i]-data2[:,i]),label=labels[i],linewidth=3)
                i = i + 1
        ax[0].legend(loc="upper right")
        ax[2].set_xlabel(r"log(Temperature)")
        ax[2].set_ylabel("Scaled / Cloudy")
	ax[0].set_title("Z = -0.3")
	ax[1].set_title("Z = 0")
	#ax[2].set_autoscale_on(False)
        #ax[2].set_aspect('equal')


	plt.savefig(outputfile)
	return data1, data2


def compare_projected_cell_emission(outputfile):
	patt = "/u/10/l/lnc2115/vega/data/Ryan/frbs/frbx"
	files = [patt+"_5kpc_HaR.cpkl",patt+"_5kpc_CIV_Scaled_ncut.cpkl",patt+"_5kpc_OVI_Scaled_ncut.cpkl",patt+"_5kpc_CIII_977_Scaled_ncut.cpkl"]
	energies = [3.028e-12,1.28e-11,1.92e-11]
	labels = ["Ha","CIV","OVI","CIII 977"]

	datain=triangle_from_frb(files,energies,labels,outputfile)
	return datain

def visualize_frb(inputfile,title,outputfile,energy=1.0,scale=False):
	if scale==True:
		frb = scale_by_energy(inputfile,energy)
	else:
		frb = cPickle.load(open(inputfile,'rb'))

	plt.imshow(np.log10(frb))
	plt.title(title)
	plt.colorbar()
	plt.savefig(outputfile)
	plt.close()
	#return

	return frb

def plot_Werk_ColDens(ax,data,ion,xkey):
	#for a given ion, plot the column density as a function of xkey
	idx = np.where((data['ion'] == ion) & (data['logNA'] > 0.1))[0]
	logNA = data['logNA'][idx]
	x = data[xkey][idx]
	l_logNA = data['l_logNA'][idx]
	e_logNA = data['e_logNA'][idx]	

	upper = np.where(l_logNA == 'u')[0]
	lower = np.where(l_logNA == 'l')[0]
	norm  = np.where(l_logNA == 'n')[0]

	ax.errorbar(x[upper],logNA[upper],yerr=0.05,uplims=True,fmt='.',color='m')
	ax.errorbar(x[lower],logNA[lower],yerr=0.05,lolims=True,fmt='.',color='m')
	ax.errorbar(x[norm],logNA[norm],yerr=e_logNA[norm],fmt='.',color='m')
	
	return
