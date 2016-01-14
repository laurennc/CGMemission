#from lauren import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import cPickle
import itertools
from werk_galaxies import galaxies

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

def plot_Werk_ColDens(ax,data,ion,xkey,xmax=500.):
	#for a given ion, plot the column density as a function of xkey
	idx = np.where((data['ion'] == ion) & (data['logNA'] > 0.1) & (data[xkey] <= xmax))[0]
	logNA = data['logNA'][idx]
	x = data[xkey][idx]
	l_logNA = data['l_logNA'][idx]
	e_logNA = data['e_logNA'][idx]	
	ID = data['ID'][idx]

##CODE FROM MUNIER
	limKwargs = {'fmt':None,'capsize':3,'elinewidth':2,'mew':0,'yerr':0.3}
	nrmKwargs = {'fmt':'s','markersize':5} 
	for i in range(len(x)):
		sSFR = galaxies[ID[i].strip()][5]/10**(galaxies[ID[i].strip()][2])
		#color = 'DodgerBlue' if sSFR > 1e-11 else 'Crimson'
		color = 'MediumBlue' if sSFR > 1e-11 else 'Crimson'
		if galaxies[ID[i].strip()][5] > 5.5:
			color = 'MediumBlue'#'LawnGreen'
		if l_logNA[i] == 'u':
			ax.errorbar(x[i],logNA[i],lolims=True,ecolor=color,**limKwargs)
		elif l_logNA[i] == 'l':
			ax.errorbar(x[i],logNA[i],uplims=True,ecolor=color,**limKwargs)
		elif l_logNA[i] == 'n':
			ax.errorbar(x[i],logNA[i],yerr=e_logNA[i],color=color,**nrmKwargs)

	#upper = np.where(l_logNA == 'u')[0]
	#lower = np.where(l_logNA == 'l')[0]
	#norm  = np.where(l_logNA == 'n')[0]

	#ax.errorbar(x[upper],logNA[upper],yerr=0.5,uplims=True,fmt=None,ecolor='m',capsize=5,elinewidth=2,mew=0)
	#ax.errorbar(x[lower],logNA[lower],yerr=0.5,lolims=True,fmt=None,ecolor='Crimson',capsize=5,elinewidth=2,mew=0)
	#ax.errorbar(x[norm],logNA[norm],yerr=e_logNA[norm],fmt='o',color='DarkOrange',ecolor='m')
	return


def plot_ion_fractions(inputfile,ax,txtDens):
	#inputfile = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/Ions/control/bertone/euvb_ion_run10_O.dat"
	#print inputfile
	data = np.genfromtxt(inputfile,skip_header=12)	
	x = data[:,0]
	i = 1
	while (i < len(data[0,:])):
		ax.plot(x,data[:,i],linewidth=1.5,label=str(i))
		i = i + 1
	ax.set_ylim(-10,1)
	#ax.text(7,-2,txtDens)
	return


def plot_ion_loop(inputPatts,outputfile,ion,xlen,ylen,special=False):
        fig,ax = plt.subplots(ylen,xlen,sharex=True,sharey=True)
        if ylen > 1:
		fig.set_size_inches(13.5,16)
	else:
		fig.set_size_inches(12,4)
        fig.subplots_adjust(hspace=0.1,wspace=0.1)
	labels = ['Factor=0.001','Factor=0.1','Factor=1','Factor=10','Factor=1000']
	j = 0
	while j < ylen:
		pattnow = inputPatts[j]
		inputfiles = [pattnow+'15_'+ion+'.dat',pattnow+'7_'+ion+'.dat',pattnow+'1_'+ion+'.dat']	
		print inputfiles
		if special == True:
			if j > 0:
				inputfiles = [pattnow+'3_'+ion+'.dat',pattnow+'2_'+ion+'.dat',pattnow+'1_'+ion+'.dat']
			#labels=["HM05","HM96","CoronalWith"]
			labels = ["All","Gals","Quasar"]
			ax[0].set_ylabel(labels[0])
			ax[1].set_ylabel(labels[1])
			ax[2].set_ylabel(labels[2])

		densLabels=['n=1','n=-3','n=-6']
		i = 0
		while i < xlen:
			if ylen > 1:
                		plot_ion_fractions(inputfiles[i],ax[j,i],densLabels[i])
                		if xlen == 0:
					ax[j,i].set_ylabel(labels[j])
			else:
				plot_ion_fractions(inputfiles[i],ax[i],densLabels[i])          
				if xlen == 0:
					ax[i].set_ylabel(labels[j])
			i = i + 1
		j = j + 1
        if ylen > 1:
		ax[0,0].set_title(densLabels[0],fontsize=18)
		ax[0,1].set_title(densLabels[1],fontsize=18)
		ax[0,2].set_title(densLabels[2],fontsize=18)
		plt.suptitle('Ion Fractions of '+ion,fontsize=22.)
		ax[0,2].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	else:
		plt.suptitle('Ion Fractions of '+ion,fontsize=22.)
		ax[2].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.savefig(outputfile)
        return

def imshow_frb(modelname,fileout):
        frbarr = np.array(cPickle.load(open(modelname,'rb')))
        #cmap = colors.ListedColormap(['Gray','HotPink','DarkTurquoise','Chartreuse'])
        bounds = [-5,1,2,3,5]
        #norm = colors.BoundaryNorm(bounds,cmap.N)
        im = plt.imshow(np.log10(frbarr),extent=(-160,160,160,-160),interpolation='none')#,vmin=-5,vmax=5,interpolation='none')#,cmap=cmap,norm=norm)
        plt.colorbar()
	plt.savefig(fileout)
	plt.close()
	#if include_colorbar==True:
        #       plt.colorbar(im)
        return


def plot_all_ions(fileout):
	data = cPickle.load(open('werk_coldens_data.cpkl','rb'))
	ions = ['HI','MgII','SiII','SiIII','SiIV','CIII','OVI']
	colors = itertools.cycle(['Indigo','b','LimeGreen','Gold','OrangeRed','MediumVioletRed','Crimson'])
	for ion in ions:
		idx = np.where((data['ion'] == ion) & (data['logNA'] > 0.1) & (data['Rperp'] <= 120) & (data['Rperp'] >= 30))[0] # & (data['l_logNA']=='n'))[0]
		upper = np.where(data['l_logNA'][idx] == 'u')[0]
		lower = np.where(data['l_logNA'][idx] == 'l')[0]
		norm = np.where(data['l_logNA'][idx] == 'n')[0]

		colornow = next(colors)
		plt.errorbar(data['Rperp'][idx][upper],data['logNA'][idx][upper],yerr=0.3,lolims=True,fmt=None,ecolor=colornow,capsize=5,elinewidth=2,mew=0,alpha=0.75)
		plt.errorbar(data['Rperp'][idx][lower],data['logNA'][idx][lower],yerr=0.3,uplims=True,fmt=None,ecolor=colornow,capsize=5,elinewidth=2,mew=0,alpha=0.75)
		plt.errorbar(data['Rperp'][idx][norm],data['logNA'][idx][norm],yerr=data['e_logNA'][idx][norm],fmt='o',ecolor=colornow,label=ion,alpha=0.75)		

	plt.xlabel('Radius (kpc)')
	plt.ylabel('Column Density')
	plt.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.,numpoints=1)
	plt.savefig(fileout,bbox_inches='tight')
	plt.close()
	return


def plot_specific_ion_fraction(inputfile,ax,ionnum,color,labelOn=False,txt='',plotDashed=False,plotDotted=False):
	data = np.genfromtxt(inputfile,skip_header=12)
	x = data[:,0]
        if labelOn:
		ax.plot(x,data[:,ionnum],linewidth=1.5,color=color,label=txt)
	elif plotDashed:
		ax.plot(x,data[:,ionnum],linewidth=1.5,linestyle='--',color=color)
	elif plotDotted:
		ax.plot(x,data[:,ionnum],linewidth=1.5,linestyle='..',color=color)
	else:
		ax.plot(x,data[:,ionnum],linewidth=1.5,color=color)

        ax.set_ylim(-10,1)
        #ax.set_ylim(0,1)
	#ax.set_xlim(4.5,7)
	return

def plot_specific_ion_EUVB_loop(inputPatt,outputfile,ion,ionnum):
	fig,ax = plt.subplots(1,3,sharey=True)
	fig.set_size_inches(12,4)
	fig.subplots_adjust(hspace=0.1,wspace=0.1)
	
	euvb = ['/g1q01/bert_ion_run','/g1q1/bert_ion_run','/g1q10/bert_ion_run']
	densPatt = ['n=1','n=-3','n=-6']
	colors = ['MediumVioletRed','Coral','Teal']

	j = 0
	while j < len(euvb):
		inputfiles=[inputPatt+euvb[j]+'15_'+ion+'.dat',inputPatt+euvb[j]+'7_'+ion+'.dat',inputPatt+euvb[j]+'1_'+ion+'.dat']
		k = 0
		while k < len(inputfiles):
			if j == (len(euvb)-1):
				plot_specific_ion_fraction(inputfiles[k],ax[j],ionnum,colors[k],labelOn=True,txt=densPatt[k])
			else:
				plot_specific_ion_fraction(inputfiles[k],ax[j],ionnum,colors[k])
			k = k + 1

		j = j + 1

	plt.suptitle('Ion Fraction of '+ion+str(ionnum),fontsize=18.)
	ax[2].legend()
	plt.savefig(outputfile)
	plt.close()	
	return

def plot_specific_ion_dens_loop(inputPatt,outputfile,ion,ionnum):
        fig,ax = plt.subplots(1,4,sharey=True)
        fig.set_size_inches(16,4)
        fig.subplots_adjust(hspace=0.1,wspace=0.1)

        euvb = ['/g1q01/bert_ion_run','/g1q1/bert_ion_run','/g1q10/bert_ion_run']
        dens = ['15_'+ion+'.dat','7_'+ion+'.dat','5_'+ion+'.dat','1_'+ion+'.dat']
	euvbPatt = ['g1q01','g1q1','g1q10']
	densPatt = ['n=1','n=-3','n=-4','n=-6']
        colors = ['MediumVioletRed','Coral','Teal']

        j = 0
        while j < len(dens):
                inputfiles=[inputPatt+euvb[0]+dens[j],inputPatt+euvb[1]+dens[j],inputPatt+euvb[2]+dens[j]]
                k = 0
                while k < len(inputfiles):
                        if j == (len(dens)-1):
				plot_specific_ion_fraction(inputfiles[k],ax[j],ionnum,colors[k],labelOn=True,txt=euvbPatt[k])
                        elif k == 2:
				plot_specific_ion_fraction(inputfiles[k],ax[j],ionnum,colors[k],plotDashed=True)
			else:
                                plot_specific_ion_fraction(inputfiles[k],ax[j],ionnum,colors[k])
                        k = k + 1
		ax[j].set_title(densPatt[j])
                j = j + 1

        plt.suptitle('Ion Fraction of '+ion+str(ionnum),fontsize=18.)
        ax[len(dens)-1].legend()
        plt.savefig(outputfile)
        plt.close()
        return

def plot_single_specific_ion(inputPatt,outputfile,ion,ionnum):
	fig,ax = plt.subplots(1,1)
	fig.set_size_inches(4,4)
	fig.subplots_adjust(hspace=0.1,wspace=0.1)
	
	ion,ionnum = 'O',6
	
	inputPatt = '/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/Ions/grid_galquas'
	euvb = ['/g1q1/bert_ion_run']
	densPatt = ['n=1','n=-3','n=-6']
	colors = ['MediumVioletRed','Coral','Teal']
	
	j = 0
	        
	inputfiles=[inputPatt+euvb[j]+'15_'+ion+'.dat',inputPatt+euvb[j]+'7_'+ion+'.dat',inputPatt+euvb[j]+'1_'+ion+'.dat']
	k = 0
	while k < len(inputfiles):
	   plot_specific_ion_fraction(inputfiles[k],ax,ionnum,colors[k],labelOn=True,txt=densPatt[k])
	                   
	   k = k + 1

	plt.suptitle('Ion Fraction of '+ion+str(ionnum),fontsize=18.)
	ax.legend()
	plt.savefig('ionfrac_Standard_O6.png')
	plt.close()
	return


def weighted_temp_dens_hist(data_sphere,weight_field,fileout):
	hden=np.log10(np.array(data_sphere["H_NumberDensity"]))
	temp=np.log10(np.array(data_sphere["Temperature"]))
	weights = np.log10(data_sphere[weight_field])
	hist = np.histogram2d(hden,temp,weights=weights,range=[[-6,2],[3, 8]])
	plt.imshow(hist[0],extent=[-6,2,3,8],interpolation='none')
	plt.colorbar()
	plt.set_cmap('jet')
	plt.savefig(fileout)
	plt.close()
	return


def plot_yt_weighted_fields(hdenfile,tempfile,fileout):
	hden = cPickle.load(open(hdenfile,'rb'))
	temp  = cPickle.load(open(tempfile,'rb'))
	plt.plot(np.log10(hden),np.log10(temp),'.',alpha=0.2,color='DodgerBlue')
	plt.xlabel('Hydrogen Density [cm^-2]')
	plt.ylabel('Temperature [K]')	
	plt.savefig(fileout)
	return


