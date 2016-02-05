import matplotlib
matplotlib.use('Agg')
import cPickle
import numpy as np
import matplotlib.pyplot as plt

model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/emis/'
model_gqs = ['g1q1']
res_keys = ['1kpc']
redshift_keys = ['z02']
znow = [0.]
ions = ['SiIV','CIII_977','CIV','OVI']
xlen,ylen = 1,4
figxlen,figylen = 3,9
fileout = 'compare_normal_hdenStepVary_z02_g1q1.png'

fig,ax = plt.subplots(ylen,xlen,sharey=True)
fig.set_size_inches(figxlen,figylen)
ax = ax.flat
i = 0

for ion in ions:
	print ion
	count = 0
	while count < len(model_gqs):
		modelnames = [model_beg+redshift_keys[count]+'/'+model_gqs[count]+'/frbx_'+res_keys[count]+'_500kpc_'+redshift_keys[count]+'_'+ion+'.cpkl',model_beg+redshift_keys[count]+'/vary_hden_step/frbx_'+res_keys[count]+'_500kpc_'+redshift_keys[count]+'_'+ion+'.cpkl']


		original = np.array(cPickle.load(open(modelnames[0],'rb'))).flat
		updated = np.array(cPickle.load(open(modelnames[1],'rb'))).flat

		ax[i].plot(np.log10(original),np.log10(updated),'b.',alpha=0.2)
		count = count + 1
	ax[i].plot(np.linspace(-8,8,8),np.linspace(-8,8,8),'k--')
	ax[i].set_title(ion)
	ax[i].set_xlim(-6,6)
	ax[i].set_ylim(-6,6)
	i = i + 1

ax[0].set_ylabel('Smaller hden steps')
ax[1].set_xlabel('Original hden steps')
		
plt.savefig(fileout)
plt.close()




		

