#from yt.mods import *
import numpy as np
import numpy
import cPickle
import matplotlib.pyplot as plt
from scipy import interpolate
from radial_data_lauren import *
import triangle

def ergs_sr_TO_raleighs(data_arr):
	return data_arr*3.304e11/79577.4715459

def LU_to_ergs(val,l_emit):
	#LU to ergs s^-1 cm^-2 arcsec^-2	
	#l_emit MUST BE IN ANGSTROMS
	return (val*4.6690e-19)/(l_emit)

def ergs_to_LU(val,l_emit):
	#ergs to LU = photons s^-1 sr^-1 cm^-2
	#l_emit MUST BE IN ANGSTROMS
	return (l_emit*val)*(2.1418e18)

def distance_from_center(x,y,z,center):
        return ((x-center[0])**2.0+(y-center[1])**2.0+(z-center[2])**2.0)**0.5

def load_Ryan_data():
	#fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"
	fn = "/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
	pf = load(fn, file_style="%s.grid.cpu%%04i") # load data
	pos = [0.40328598,0.47176743,0.46131516]
	rad = 108.0/pf['kpc']
	data = pf.h.sphere(pos,rad)
	return pf, data

def scale_by_metallicity(values,assumed_Z,wanted_Z):
	wanted_ratio = (10.**(wanted_Z))/(10.**(assumed_Z))
	return values*wanted_ratio

def scale_by_energy(frbfile,energy):
	data = cPickle.load(open(frbfile,'rb'))
        return data*(5.7e-18)*(1./1.87e-12)/(4.*np.pi*energy)

def make_Cloudy_table(table_index):
	#hden_n_bins, hden_min, hden_max = 35, -6, 1
	#T_n_bins, T_min, T_max = 151, 3, 8
	#hden_n_bins, hden_min, hden_max = 15, -6, 1
	#hden_n_bins, hden_min, hden_max = 81, -6, 2
	hden_n_bins, hden_min, hden_max = 17, -6, 2
	T_n_bins, T_min, T_max = 51, 3, 8
	patt = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/grid_galquas/emis/z02/g1q1/g1q1_run%i.dat"
	#patt = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/test_hden_step/g1q1_run%i.dat"

	hden=numpy.linspace(hden_min,hden_max,hden_n_bins)
	T=numpy.linspace(T_min,T_max, T_n_bins)
	table = np.zeros((hden_n_bins,T_n_bins))
	for i in range(hden_n_bins):
		table[i,:]=[float(l.split()[table_index]) for l in open(patt%(i+1)) if l[0] != "#"]
	return hden,T,table

def make_ion_table(ion,number):
	#hden_n_bins,hden_min,hden_max = 15, -6, 1
	hden_n_bins, hden_min, hden_max = 17, -6, 2
	T_n_bins, T_min, T_max = 51, 3, 8
	#was 10 before HI run
	#patt = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/Ions/grid_galquas/g1q01/bert_ion_run%i_"+ion+".dat"
	patt = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/Ions/grid_galquas/z02/g1q1/bert_ion_run%i_"+ion+".dat"
	hden=numpy.linspace(hden_min,hden_max,hden_n_bins)
        T=numpy.linspace(T_min,T_max, T_n_bins)
        table = np.zeros((hden_n_bins,T_n_bins))
        for i in range(hden_n_bins):
                table[i,:]=[float(l.split()[number]) for l in open(patt%(i+1)) if l[0] != "#"]
        return hden,T,table
 
def make_SB_profile(filex,filey,filez,xL,z=0.):
        #xL = np.arange(-20,20)*10.0
        #xL = np.arange(-200,200,0.5)
	#xL = np.arange(-100,100,0.25)
	#xL = np.linspace(-160,160,320)#[0:320]
	xL, yL = np.meshgrid(xL,xL)
        #r = abs(xL+1j*yL)

        frbx = cPickle.load(open(filex,'rb'))
	frby = cPickle.load(open(filey,'rb'))
	frbz = cPickle.load(open(filez,'rb'))

        rp_x = radial_data(frbx/(1.+z)**4.0,x=xL,y=yL)
        rp_y = radial_data(frby/(1.+z)**4.0,x=xL,y=yL)
        rp_z = radial_data(frbz/(1.+z)**4.0,x=xL,y=yL)

        rp_mean = (rp_x.mean + rp_y.mean + rp_z.mean)/3.0
        rp_median  = (rp_x.median + rp_y.median + rp_z.median)/3.0
	rp_max = map(max,rp_x.max,rp_y.max,rp_z.max)
	rp_min = map(min,rp_x.min,rp_y.min,rp_z.min)
	rp_std = (rp_x.std + rp_y.std + rp_z.std)/3.0
        return rp_x.r, rp_mean, rp_median, rp_max, rp_min, rp_std



def make_SB_profile_OLD(filex,filey,filez,energy):
	xL = np.arange(-20,20)*10.0
	xL, yL = np.meshgrid(xL,xL)
	r = abs(xL+1j*yL)

	frbx = cPickle.load(open(filex,'rb'))
	frby = cPickle.load(open(filey,'rb'))
	frbz = cPickle.load(open(filez,'rb'))
	
	frbx = frbx*(5.7e-18)*(1./1.87e-12)/(4.*np.pi*energy)
	frby = frby*(5.7e-18)*(1./1.87e-12)/(4.*np.pi*energy)
	frbz = frbz*(5.7e-18)*(1./1.87e-12)/(4.*np.pi*energy)

	rp_Ralx = radial_data(frbx,x=xL,y=yL)
	rp_Raly = radial_data(frby,x=xL,y=yL)
	rp_Ralz = radial_data(frbz,x=xL,y=yL)
	
	rp_mean = (rp_Ralx.mean + rp_Raly.mean + rp_Ralz.mean)/3.0
	rp_median  = (rp_Ralx.median + rp_Raly.median + rp_Ralz.median)/3.0

	return rp_Ralx.r, rp_mean, rp_median

def triangle_from_frb(files,energies,labels,outputfile):
	frb = scale_by_energy(files[0],energies[0])
	datain = np.zeros((len(frb.flatten()),len(files)))
	datain[:,0] = np.log10(frb.flatten())
	i = 1
	while i < len(files):
		#Right now, I made CIII projections properly so they don't need to be scaled....
		if i < 3:
			frb = scale_by_energy(files[i],energies[i])
		else:
			frb = cPickle.load(open(files[i],'rb'))
		idx = np.where(frb.flatten() == 0.0)
		datain[:,i] = np.log10(frb.flatten())
		if len(idx[0]) > 0:
			datain[idx,i] = 12.0
		i = i + 1
	#print datain

	#triangle.corner(datain,labels=labels,range=(-1,6)).savefig(outputfile)
	triangle.corner(datain,labels=labels).savefig(outputfile)
	return datain

def emission_interpolation(field_idx,H_N,Temperature):
	from scipy import interpolate
	hden,temp,table = make_Cloudy_table(field_idx)
	hden,temp = np.meshgrid(hden,temp)
	pts = np.array((hden.ravel(),temp.ravel())).T
	sr = table.T.ravel()
	bl = interpolate.LinearNDInterpolator(pts,sr)
	dia = bl(H_N,Temperature)
	idx = np.isnan(dia)
	dia[idx] = -200.0
	return dia

def plot_scatter_percentile(ax1,data,x,y,percentile,symbol,maxr=500.):
        #data is the frb array data
        x,y = np.meshgrid(x,y)
        working_mask = np.ones(data.shape,bool)
        r = abs(x+1j*y)
        rmax = r[working_mask].max()
        dr = np.abs([x[0,0] - x[0,1]]) #* annulus_width
        radial = np.arange(maxr/dr)*dr + dr/2.
        nrad = len(radial)
        #mslope = (0.001-0.007)/nrad
        mslope = (0.0075-0.07)/nrad

        for irad in range(nrad):
                minrad = irad*dr
                maxrad = minrad + dr
                thisindex = (r>=minrad) * (r<maxrad) * working_mask
                alphahere = mslope*irad + 0.07
                ax1.plot(r[thisindex],np.log10(data[thisindex]),symbol,alpha=alphahere)
                rhere = r[thisindex].flatten()
                datahere = data[thisindex].flatten()
                meanhere = data[thisindex].mean()
                lenperc = int(len(datahere)*percentile)
                idSort = np.argsort(datahere)[::-1]
                wanted = idSort[0:lenperc]
                if len(wanted) == 0:
                        wanted = np.where(datahere == datahere.max())
                ax1. plot(rhere[wanted],np.log10(datahere[wanted]),'g.',alpha=0.2)
                #ax1.plot(rhere[wanted],np.log10(datahere[wanted]),symbol,markersize=3.5)#,markersize=0.75)                             
        return


def logU_to_hden(logU):
	return -1.*logU-5.9694002780340494


def find_covering_fraction(frb_filename,tempname,SB_lims,znow):
	frb = np.array(cPickle.load(open(frb_filename,'rb')))
	temperature = np.array(cPickle.load(open(tempname,'rb')))
	
	idx = np.where(np.log10(temperature) < 5.2)
	test = np.where((idx[1] > 130) & (idx[0] < 170) & (idx[1] < 170) & (idx[0] > 130))
	actual = [idx[0][test],idx[1][test]]
	#print 'Min value is: '+str(np.log10((frb[actual].flatten()/(1.+znow)**4.0)).min())
	#print 'Avg value is: '+str(np.log10(np.average(frb[actual].flatten()/(1.+znow)**4.0)))
	frb[actual] = -50
	testing = len(np.where(frb < -49.)[0].flatten())/float(len(frb.flatten()))
        #print 'Fraction of pixels in the disk: '+str(testing)

	frb = frb.flatten()
	frb = np.log10(frb/((1.+znow)**4.0))
	
	fractions = []
	for lim in SB_lims:
		idx = np.where(frb > lim)[0]
		frac_temp = float(len(idx))/float(len(frb))
		fractions = np.append(fractions,frac_temp)
	return fractions

def make_radius_array():
        xL = np.linspace(-160,160,320)
        maxnow = 160.
        xL, yL = np.meshgrid(xL,xL)
        r = abs(xL+1j*yL)
        dr = np.abs([xL[0,0] - xL[0,1]])
        radial = np.arange(maxnow/dr)*dr + dr/2
        nrad = len(radial)
        return r, dr, nrad


