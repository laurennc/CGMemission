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
ions = ['HI','MgII','SiII','SiIII','SiIV','CIII','OVI']

print 'g = 0.01'
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g01q01/frb"
models = np.append(models,Model(ions,patt,model_mid,0.01,0.01,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g01q1/frb"
models = np.append(models,Model(ions,patt,model_mid,0.01,1.,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g01q10/frb"
models = np.append(models,Model(ions,patt,model_mid,0.01,10.,l,u,model_width,nbins,ndraws))#k))
#patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g01q0.05/frb"
#models = np.append(models,Model(ions,patt,model_mid,0.01,0.05,l,u,model_width))#k))
#patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g01q0.5/frb"
#models = np.append(models,Model(ions,patt,model_mid,0.01,0.5,l,u,model_width))#k))
#patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g01q2/frb"
#models = np.append(models,Model(ions,patt,model_mid,0.01,2.0,l,u,model_width))#k))

print 'g = 1'
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g1q01/frb"
models = np.append(models,Model(ions,patt,model_mid,1.,0.01,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g1q1/frb"
models = np.append(models,Model(ions,patt,model_mid,1.,1.,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g1q10/frb"
models = np.append(models,Model(ions,patt,model_mid,1.,10.,l,u,model_width,nbins,ndraws))#k))
#patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g1q0.05/frb"
#models = np.append(models,Model(ions,patt,model_mid,1.,0.05,l,u,model_width))#k))
#patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g1q0.5/frb"
#models = np.append(models,Model(ions,patt,model_mid,1.,0.5,l,u,model_width))#k))
#patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g1q2/frb"
#models = np.append(models,Model(ions,patt,model_mid,1.,2.,l,u,model_width))#k))

print 'g = 10'
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g10q01/frb"
models = np.append(models,Model(ions,patt,model_mid,10.,0.1,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g10q1/frb"
models = np.append(models,Model(ions,patt,model_mid,10.,1.,l,u,model_width,nbins,ndraws))#k))
patt = "/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/coldens/grid_galquas/g10q10/frb"
models = np.append(models,Model(ions,patt,model_mid,10.,10.,l,u,model_width,nbins,ndraws))

print 'starting plots'
count = 0
iax = 130
fig = plt.figure(figsize = (12,4))
markers = itertools.cycle(['^','o','s'])#,'v','p','D'])
colors = itertools.cycle(['b','g','r'])#,'c','m','k'])
for count in range(len(models)):
        print count
	if count %3 == 0.0:
                iax = iax+1
                ax1 = fig.add_subplot(iax)
                ax1.set_title('g='+str(models[count].g))
                ax1.set_xticks(np.arange(0, len(ions), 1.0))
                ax1.set_xticklabels(ions)
        markernow= next(markers)
        colornow = next(colors)
        for idl in range(len(ions)):
                ax1.plot(idl,models[count].lines[idl].total_likelihood,color=colornow,marker=markernow,label='q='+str(models[count].q),alpha=0.6)

        ax1.set_xlim(-0.2,len(ions)+0.02)
plt.savefig('model_likelihoods_nbins500_ndraws10_detectonly.png')
plt.close()





###################THIS WILL ALL BE MORE USEFUL WHEN I HAVE MORE MODELS########################################33
#iax = 331
#fig = plt.figure(figsize=(13.5,13.5))

#for idl in range(len(ions)):
#	ax1 = fig.add_subplot(iax)
#	print ions[idl]
#	for idx in range(len(models)):
#		#print models[idx].g,models[idx].q,models[idx].lines[idl].ion,models[idx].lines[idl].likelihood
#		sup = ax1.scatter(models[idx].g,models[idx].q,c=models[idx].lines[idl].likelihood,marker='o')#,markersize=8.)
#		ax1.set_title(ions[idl])
#	iax = iax+ 1
#
#fig.subplots_adjust(right=0.75)
#cbar_ax = fig.add_axes([0.8, 0.15, 0.05, 0.7])
#fig.colorbar(sup, cax=cbar_ax,ticks=[-200, -100, -50,0])
#
#plt.gcf().subplots_adjust(bottom=0.15)
#plt.savefig('model_fits.png')
#plt.close()



