#Class designed to hold the Lines of each model for plotting and other calculations
#Created by: Lauren
#Created on: 10/27/2014
import matplotlib
matplotlib.use('Agg')
from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import itertools
from Line_ORG import *
from Model_ORG import *

qs = [0.01,1,10.]
gs = [0.01,1,10.]
l,u,k = 10.,20.,8.
model_width = 0.0#0.5
nbins, ndraws = 500,10
models = []

model_mid = '_1kpc_z02_'
#ions = ['CIII','CIV','OVI','MgII','SiII','SiIII','SiIV','HI']
#ions = ['CIII','OVI','MgII','SiII','SiIII','SiIV','HI']
ions = ['HI','MgII','SiII','SiIII','SiIV','CIII']#,'OVI']

print 'g = 0.01'
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g01q0.05/frb"
models = np.append(models,Model(ions,patt,model_mid,0.01,0.05,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g01q01/frb"
models = np.append(models,Model(ions,patt,model_mid,0.01,0.01,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g01q0.1/frb"
models = np.append(models,Model(ions,patt,model_mid,0.01,0.1,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g01q0.5/frb"
models = np.append(models,Model(ions,patt,model_mid,0.01,0.5,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g01q1/frb"
models = np.append(models,Model(ions,patt,model_mid,0.01,1.,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g01q2/frb"
models = np.append(models,Model(ions,patt,model_mid,0.01,2.0,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g01q10/frb"
models = np.append(models,Model(ions,patt,model_mid,0.01,10.,l,u,model_width,nbins,ndraws))#k))


print 'g = 1'
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g1q0.05/frb"
models = np.append(models,Model(ions,patt,model_mid,1.,0.05,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g1q01/frb"
models = np.append(models,Model(ions,patt,model_mid,1.,0.01,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g1q0.1/frb"
models = np.append(models,Model(ions,patt,model_mid,1.,0.1,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g1q0.5/frb"
models = np.append(models,Model(ions,patt,model_mid,1.,0.5,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g1q1/frb"
models = np.append(models,Model(ions,patt,model_mid,1.,1.,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g1q2/frb"
models = np.append(models,Model(ions,patt,model_mid,1.,2.,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g1q10/frb"
models = np.append(models,Model(ions,patt,model_mid,1.,10.,l,u,model_width,nbins,ndraws))#k))

print 'g = 10'
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g10q01/frb"
models = np.append(models,Model(ions,patt,model_mid,10.,0.1,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g10q1/frb"
models = np.append(models,Model(ions,patt,model_mid,10.,1.,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g10q10/frb"
models = np.append(models,Model(ions,patt,model_mid,10.,10.,l,u,model_width,nbins,ndraws))



likelihoods  = [model.model_probability for model in models]
x = np.arange(len(models)) + 0.5
ticklabels = ['g01q0.05','g01q01','g01q0.1','g01q0.5','g01q1','g01q2','g01q10','g1q0.05','g1q01','g1q0.1','g1q0.5','g1q1','g1q2','g1q10','g10q01','g10q1','g10q10']

plt.plot(x,likelihoods,'s',color='DeepPink')
plt.xticks(x, ticklabels, rotation=45)
plt.tight_layout()
plt.savefig('total_model_evaluation_sinOVI_dectectonly.png')
plt.close()


