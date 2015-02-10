import matplotlib
matplotlib.use('Agg')
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import itertools

data = cPickle.load(open('werk_coldens_data.cpkl','rb'))
ions = ['HI','MgII','SiII','SiIII','SiIV','CIII','OVI']
ion_x = {'HI':0.5,'MgII':1.5,'SiII':2.5,'SiIII':3.5,'SiIV':4.5,'CIII':5.5,'OVI':6.5}

count = 0
plot_count = 0 

fig,axs = plt.subplots(4,4)
fig.set_size_inches(13.5,13.5)
ax2 = axs.flat

all_IDs = np.unique(data['ID'])

while count < len(all_IDs):
	iax = 0
	fig,axs = plt.subplots(4,4)
	fig.set_size_inches(13.5,13.5)
	ax2 = axs.flat	
	while iax < len(ax2):
		ID = all_IDs[count]
		idx = np.where((data['ID'] == ID) & (data['logNA'] > 0.1) & (data['Rperp'] <= 120) & (data['Rperp'] >= 30))[0] # & (data['l_logNA']=='n'))[0]
		
		ions_here = data['ion'][idx]
		coldens = data['logNA'][idx]
		data_label = data['l_logNA'][idx]
		err_here = data['e_logNA'][idx]
		
		xarr,yarr,yerr,ylabel,j = [],[],[],[],0
		while j<len(ions_here):
			if ions_here[j] in ions:
				xarr = np.append(xarr,ion_x[ions_here[j]])
				yarr = np.append(yarr,coldens[j])
				yerr = np.append(yerr,err_here[j])
				ylabel = np.append(ylabel,data_label[j])
			j = j + 1

		upper = np.where(ylabel == 'u')[0]
       		lower = np.where(ylabel == 'l')[0]
        	norm  = np.where(ylabel == 'n')[0]
		
		if len(xarr) == 0:
			iax = iax -1
		else:
			if len(upper)>0:
				ax2[iax].errorbar(xarr[upper],yarr[upper],yerr=0.3,uplims=True,fmt=None,ecolor='m',capsize=5,elinewidth=2,mew=0)
			if len(lower)>0:
				ax2[iax].errorbar(xarr[lower],yarr[lower],yerr=0.3,lolims=True,fmt=None,ecolor='Crimson',capsize=5,elinewidth=2,mew=0)
			if len(norm)>0:
				ax2[iax].errorbar(xarr[norm],yarr[norm],yerr=yerr[norm],fmt='o',color='DarkOrange',ecolor='m')
			
			ax2[iax].set_xticks(np.arange(0.5,len(ions)+0.5,1.0))
			ax2[iax].set_xticklabels(ions)
			ax2[iax].set_title(ID)
			#ax2[iax].set_ylabel('log(N_ion)')
		
		iax = iax + 1
		count = count + 1
		if count == len(all_IDs):
			break

	plt.subplots_adjust(hspace=0.4)	
	plt.savefig('iondata_by_galaxy'+str(plot_count)+'.png')
	plt.close()

	plot_count = plot_count + 1


