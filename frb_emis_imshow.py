#Code to plot the frbs and radial profiles or my emission predictions
#Created by: Lauren
#Created on: 12/11/2014
import matplotlib
matplotlib.use('Agg')
import cPickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import itertools
import matplotlib.gridspec as gridspec
import re
from lauren_holding import make_SB_profile
from radial_data_lauren import *
from plotting_routines import *
from matplotlib import colors
import brewer2mpl as brew

############FUNCTIONS FOR THE PLOTTING####################################

def basic_profile_vals(modelnames,xL,z=0.):
        rpr,rpmean,rpmedian,rpmax,rpmin,rp_std = make_SB_profile(modelnames[0],modelnames[1],modelnames[2],xL,z=z)
        return rpr,rpmean,rpmedian,rpmax,rpmin

def make_radius_array(key):
	if key=='1kpc':
        	xL = np.linspace(-160,160,320)
        elif key=='5kpc':
		xL = np.linspace(-160,160,64)
	elif key=='25kpc':
		xL = np.linspace(-160,160,13)
	else:
		print 'Failure in the key'
	maxnow = 160.
        xL2, yL = np.meshgrid(xL,xL)
        r = abs(xL2+1j*yL)
        dr = np.abs([xL2[0,0] - xL2[0,1]])
        radial = np.arange(maxnow/dr)*dr + dr/2
        nrad = len(radial)
        return xL,r, dr, nrad

def plot_scatter_points_percentile(ax,frbarr,r,dr,nrad,percentile,key,z=0.,comoving=False):
        mslope = (0.0075-0.07)/nrad
        for irad in range(nrad):
		minrad = irad*dr
        	maxrad = minrad + dr
        	thisindex = (r>=minrad) * (r<maxrad)
        	alphahere = mslope*irad + 0.07
		rnow,frbnow = r[thisindex],np.log10(frbarr[thisindex]/(1.+z)**4.0)	
	
		if comoving:
			rnow = rnow*(1.+z)

		lims = [10,3,2,1,-10]
                colors = ['Chartreuse','DarkTurquoise','HotPink','Grey']
                i = 0
		while i < len(lims)-1:
			wanted = np.where((frbnow < lims[i]) & (frbnow >= lims[i+1]))[0]
			if key=='1kpc':
				ax.plot(rnow[wanted],frbnow[wanted],'.',color=colors[i],alpha=alphahere)
			else:
				ax.plot(rnow[wanted],frbnow[wanted],'.',color=colors[i],alpha=0.37)
			i = i + 1

		#rhere = r[thisindex].flatten()
                #datahere = frbarr[thisindex].flatten()
                #meanhere = frbarr[thisindex].mean()
                #lenperc = int(len(datahere)*percentile)
                #idSort = np.argsort(datahere)[::-1]
                #wanted = idSort[0:lenperc]
                #if len(wanted) == 0:
                #        wanted = np.where(datahere == datahere.max())
                #ax.plot(rhere[wanted],np.log10(datahere[wanted]),'g.',alpha=0.5)
        return

def full_scatter_plot(modelnames,ion,ax,res_key,max_r,percentile,znow,comoving=False):
	xL,r,dr,nrad = make_radius_array(res_key)
	##EDITED HERE FROM 0 TO 1
	frbarr = np.array(cPickle.load(open(modelnames[0],'rb')))
	plot_scatter_points_percentile(ax,frbarr,r,dr,nrad,percentile,res_key,z=znow,comoving=comoving)
	rpr,rpmean,rpmedian,rpmax,rpmin = basic_profile_vals(modelnames,xL,z=znow)
	
	idr = np.where(rpr <= 160.0)
	if comoving:
		ax.plot(rpr[idr]*(1.+znow),np.log10(rpmedian[idr]),color='k',linewidth=2.2)
	else:
		ax.plot(rpr[idr],np.log10(rpmedian[idr]),color='k',linewidth=2.2)#previously Goldenrod	
	#print np.max(rpr[idr]),np.max(np.log10(rpmedian[idr]))
	ax.set_xticks(range(0,160,30))
	#ax.set_xlim(0,160)
	ax.axis([0.,160.,-6.,6.5])
	return

def plot_frb(modelname,ax,z=0.,include_colorbar=False,obs_colors=True):
	frbarr = np.array(cPickle.load(open(modelname,'rb')))
	if obs_colors:
		cmap = colors.ListedColormap(['Gray','HotPink','DarkTurquoise','Chartreuse'])
		bounds = [-5,1,2,3,5]	
		norm = colors.BoundaryNorm(bounds,cmap.N)
		im = ax.imshow(np.log10(frbarr/(1.+z)**4.0),extent=(-160,160,160,-160),vmin=-5,vmax=5,interpolation='none',cmap=cmap,norm=norm,origin='lower')
	else:
		bmap = brew.get_map('PuRd','Sequential',9)
		#bmap = brew.get_map('PRGn','Diverging',9)
		cmap = bmap.get_mpl_colormap(N=1000, gamma=2.0)
		im = ax.imshow(np.log10(frbarr/(1.+z)**4.0),extent=(-160,160,160,-160),vmin=-5,vmax=5,interpolation='none',cmap=cmap,origin='lower')

	#if include_colorbar==True:
	#	plt.colorbar(im)
	return im

def add_HI_contour(ax,HIfrb,ncontours):
	xL,r,dr,nrad = make_radius_array('1kpc')
	X,Y = np.meshgrid(xL,xL)
	CS = ax.contour(X,Y,HIfrb,[15,18],color='k')
	#plt.clabel(CS, fontsize=9, inline=1)
	return

def werk_cloudy_ions(emis,gals,ion):
        blues = np.array([0.,0.])
        reds = np.array([0.,0.])

        i,j = 0,0

        ##set it up for subplots since I'll be making a nine panel plot

        while i < len(emis['ID']):
                gal_id = emis['ID'][i]
                #print gal_id
                #gal_id = gal_id[0:10]
                idx = np.where(gals['ID'] == gal_id)[0]
                sSFR = gals['SFR'][idx]/np.power(10.,gals['mStar'][idx])
                color = 'DodgerBlue' if sSFR > 1e-11 else 'Crimson'
                if color == 'DodgerBlue':
                        temp = np.array([0.,0.])
                        temp[0],temp[1] = gals['Rperp'][idx],emis[ion+'min'][j]
                        blues = np.vstack((blues,temp))
                        temp[1] = emis[ion+'max'][j]
                        blues = np.vstack((blues,temp))
                if color == 'Crimson':
                        temp = np.array([0.,0.])
                        temp[0],temp[1] = gals['Rperp'][idx],emis[ion+'min'][j]
                        reds = np.vstack((reds,temp))
                        temp[1] = emis[ion+'max'][j]
                        reds = np.vstack((reds,temp))

                i = i + 4
                j = j + 1


        blues = np.delete(blues,0,0)
        reds = np.delete(reds,0,0)
        return blues,reds

def plot_connected_minmax(ax,blues,reds):
        i = 0
        while i < len(blues):
                ax.plot(blues[i:i+2,0],blues[i:i+2,1],'.',color='DodgerBlue')
                ax.plot(blues[i:i+2,0],blues[i:i+2,1],color='DodgerBlue',lw=2.0)
                i = i + 2

        i = 0
        while i < len(reds):
                ax.plot(reds[i:i+2,0],reds[i:i+2,1],'.',color='Crimson')
                ax.plot(reds[i:i+2,0],reds[i:i+2,1],color='Crimson',lw=2.0)
                i = i + 2

        return



#########################################################################

#ions = ['SiIV','CIV','OVI']
#ions = ['CIII_977','CIV','OVI']

no_axes_labels = False
obs_colors = False

model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/final/emis/' ##CHANGED FORM Z02
#model_beg = '/u/10/l/lnc2115/vega/repos/CGMemission/bertone_frbs/emis/grid_galquas/'

##ASTROFEST PARAMETERS
#model_gqs = ['g1q1']
#res_keys = ['1kpc']
#redshift_keys = ['z0']
#znow = [0.]
#ions = ['SiIV','CIII_977','CIV','OVI']


##REDSHIFT EVOLUTION PARAMETERS
model_gqs = ['g1q1','g1q1','g1q1']
res_keys = ['1kpc','1kpc','1kpc']
redshift_keys = ['z0','z05','z1']
#znow = [0.,0.5,1.0]
znow = [0.,0.,0.]
ions = ['CIII_977','CIV','OVI']	

max_r = 160.
percentile = 0.01
ncontours = 4

fileout = 'frb_emis_theory_nozscaling.png'

xlen,ylen = 3,3
fig,ax = plt.subplots(ylen,xlen,sharey=True)
fig.set_size_inches(12,12)
#fig.set_size_inches(24,6)
#gs1 = gridspec.GridSpec(1, 4)
#gs1.update(wspace=0.025, hspace=0.05)
ax = ax.flat
i = 0


for ion in ions:
	print ion
	count = 0
	while count < len(model_gqs):
		print ion, count, model_gqs[count],i
		modelnames = [model_beg+redshift_keys[count]+'/'+model_gqs[count]+'/frbx_'+res_keys[count]+'_500kpc_'+redshift_keys[count]+'_'+ion+'.cpkl'] 

		im_out = plot_frb(modelnames[0],ax[i],z=znow[count],obs_colors=obs_colors)

		if no_axes_labels:
			ax[i].set_title(model_gqs[count])
			plt.axis('on')
			ax[i].set_xticklabels([])
			ax[i].set_yticklabels([])
			ax[i].set_yticklabels([])
			ax[i].set_adjustable('box-forced')		

		#ax[i].set_title(ion)		
		i = i + 1
		count = count + 1	

fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.9, 0.1, 0.05, 0.5])
fig.colorbar(im_out, cax=cbar_ax)

for j in range(len(model_gqs)):
        #ax[j].set_title(model_gqs[j])
	ax[j].set_title(redshift_keys[j])
	
for k in range(len(ions)):
        ax[k*xlen+2].text(75,-100,ions[k])

#plt.tight_layout()
plt.savefig(fileout)#,transparent=True)
plt.close()



