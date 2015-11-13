import matplotlib
matplotlib.use('Agg')
import cPickle
import numpy as np
import matplotlib.pyplot as plt
field = 'Temperature'
weight = 'Density'
filein = 'bertone_frbs/final/basic/z1/g1q1/'+'frbx_1kpc_500kpc_z1_'+field+'_'+weight+'.cpkl'
data = cPickle.load(open(filein,'rb'))
idx = np.where(np.log10(data) < 5.2)

model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/emis/'
model_gqs = ['g1q1']
res_keys = ['1kpc']
redshift_keys = ['z1']
znow = [0.]
ions = ['SiIV','CIII_977','CIV','OVI']
ion = 'OVI'
count = 0
modelnames = [model_beg+redshift_keys[count]+'/'+model_gqs[count]+'/frbx_'+res_keys[count]+'_500kpc_'+redshift_keys[count]+'_'+ion+'.cpkl']
frbarr = np.log10(np.array(cPickle.load(open(modelnames[0],'rb'))))
#frbarr[idx] = -50.

test = np.where((idx[1] > 130) & (idx[0] < 170) & (idx[1] < 170) & (idx[0] > 130))
actual = [idx[0][test],idx[1][test]]
frbarr[actual] = -50

im = plt.imshow(frbarr,extent=(-160,160,160,-160),vmin=-5,vmax=5,interpolation='none',origin='lower')
plt.colorbar()
plt.savefig('finding_disk_'+redshift_keys[count]+'.png')
plt.close()



