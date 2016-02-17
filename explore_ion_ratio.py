##ok so what are the different things I want to do?
## plot ion ratios by density
## plot ion ratios by temperature
## perhaps a radial plot where temperature indicated
##      --- but I'm not sure how that would work with opacity...
## phase diagram colored by ion ratio
import matplotlib
matplotlib.use('Agg')
import cPickle
import numpy as np
import matplotlib.pyplot as plt

def load_emis_file(ion,model_gq):
	model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/emis/z02/'
	model_mid = '/frbx_1kpc_500kpc_z02_'
	model = model_beg+model_gq+model_mid+ion+'.cpkl'
	return np.array(cPickle.load(open(model,'rb')))

def load_coldens_file(ion,model_gq):
        model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/coldens/z02/'
        model_mid = '/frbx_1kpc_500kpc_z02_'
        model = model_beg+model_gq+model_mid+ion+'dens.cpkl'
        return np.array(cPickle.load(open(model,'rb')))

def load_basic_file(field):
	model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/basic/z02/g1q1/'
        model_mid = '/frbx_1kpc_500kpc_z02_'
        model = model_beg+model_mid+field+'.cpkl'
        return np.array(cPickle.load(open(model,'rb')))

def basic_plot(ax,xfield,yfield1,yfield2):
	ion_frac = np.log10(yfield1/yfield2)
	ax.plot(np.log10(xfield),ion_frac,'k.',alpha=0.01)
	ax.set_xticks(np.arange(round(np.log10(xfield).min(),1), round(np.log10(xfield).max()+1,1), 1.0))
	return

def basic_imshow(ax,xfield,yfield1,yfield2):
	ion_frac = yfield1/yfield2
	ion_frac = np.log10(ion_frac.flat)
	#hist,xedges,yedges = np.histogram2d(np.log10(xfield.flat),ion_frac,range=[[4,6.5],[-3,3]],bins=100)
	#hist,xedges,yedges = np.histogram2d(np.log10(xfield.flat),ion_frac,range=[[-5,-1],[-4,2]],bins=100)
	hist,xedges,yedges = np.histogram2d(np.log10(xfield.flat),ion_frac,range=[[18,22],[-6,3]],bins=100)

	extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
	#ax.imshow(np.log10(hist.T),origin='lower',extent=[4,6.5,-3,3],interpolation='nearest',cmap='jet',vmin=-3.,vmax=3.) 
	#ax.imshow(np.log10(hist.T),origin='lower',extent=[-5,-1,-4,2],interpolation='nearest',cmap='jet',vmin=-3.,vmax=3.)
	ax.imshow(np.log10(hist.T),origin='lower',extent=[18,22,-6,3],interpolation='nearest',cmap='jet',vmin=-3.,vmax=3.)
	#ax.set_xlim(4,6.5)
	ax.set_aspect('auto')
	return


xlen,ylen = 3,2 #6#2,6
fig,ax = plt.subplots(ylen,xlen,sharex=True,sharey=True)
fig.set_size_inches(9,6)
ax = ax.flat

fileout = 'NH_CIII_OVI_imshow.png'
#xfield = 'Temperature_Density'
xfield = 'H_NumberDensity'
#xfield = 'Density'
basic = load_basic_file(xfield)

model_gqs = ['g1q01','g1q1','g1q10']
ion1,ion2,ion3 = 'CIII_977','OVI','CIII'
#ion1,ion2,ion3 = 'SiIII','SiIII','SiII'
count = 0

for i in range(len(ax)):
	print i,count
	if i < 3:
		ion_top = load_emis_file(ion1,model_gqs[count])
		ion_btm = load_emis_file(ion2,model_gqs[count])
		#basic_plot(ax[i],basic,ion_top,ion_btm)
		basic_imshow(ax[i],basic,ion_top,ion_btm)
		ax[i].set_title(model_gqs[i])
	else:
		ion_top = load_coldens_file(ion3,model_gqs[count])
		ion_btm = load_coldens_file(ion2,model_gqs[count])
		#basic_plot(ax[i],basic,ion_top,ion_btm)
		basic_imshow(ax[i],basic,ion_top,ion_btm)		

	count = count + 1
	if count == 3:
		count = 0


ax[0].set_ylabel('Emission '+ion1+'/'+ion2)
ax[3].set_ylabel('Absorption '+ion3+'/'+ion2)
ax[4].set_xlabel('N_H')
#ax[4].set_xlabel('Temperature')
#ax[4].set_xlabel('Density')
plt.savefig(fileout)
plt.close()




