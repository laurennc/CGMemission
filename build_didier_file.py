import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from yt.mods import *
import numpy as np
import sys
from astropy.cosmology import WMAP5 as cosmo
from astropy import units as u
from astropy.cosmology import z_at_value
from scipy.interpolate import interp1d
from astropy.io import fits

def set_central_pixel(ds,pos):
	idx = np.where((ds['x']==pos[0])&(ds['y']==pos[1])&(ds['z']==pos[2]))[0]
	return idx

def build_region(center,width,thickness,axis):
	axis = fix_axis(axis)
	
	center = np.array(center)
	LE,RE = center.copy(),center.copy()
	LE[axis] -= thickness/2.0
	RE[axis] += thickness/2.0
	
	area_axes = [0,1,2]
	i = area_axes.index(axis)
	del area_axes[i]
	
	LE[area_axes] -= width/2.0
	RE[area_axes] += width/2.0
	
	return pf.h.region(center, LE, RE)

def build_interpolator(filein):
	data = np.genfromtxt(filein)
	return interp1d(data[:,0],data[:,1])

 

#fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
interp_file = "redshift_lumindist.dat"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

val, pos = pf.h.find_max('Density')
rad = 320./pf['kpc']
thickness = 500./pf['kpc']
axis = ['z']

region = build_region(pos,rad,thickness,axis)
z_cos = 1.0
#los --> choosing the x-axis to be my los
los = (0,0,1)
#ra = (1,0,0)
#dec = (0,1,0)


#################################################
### BUILDING THE FIELDS THAT I NEED #############
### FIELDS ARE LISTED BELOW #####################
#################################################

##fields to be exported
x,y,z = [],[],[]
ra,dec = [],[]
z_geo,z_app = [],[]
dlum = []
vlos = []
T,rho,Z = [],[],[]
nHI,ne = [],[]
P,dx = [],[]
emiss_HI,emiss_OVI,emiss_CIV = [],[],[]

##let's projection along the z direction with pos as the reference value
##NOTE: I COULD MAKE THIS MORE GENERIC BUT FOR EASE NOW LET'S DO IT THIS WAY
x = (pos[0]-region['x'])*pf['Mpccm']
y = (pos[1]-region['y'])*pf['Mpccm']
z = (pos[2]-region['z'])*pf['Mpccm']

#easy basic quantities
vlos = region['z-velocity']*1.e-5 ##km/s

T,P,rho = region['Temperature'], region['Pressure'], region['H_NumberDensity']
#T = [K], P = [dyne/cm^2], rho = [atoms/cm^3]

emiss_HI,emiss_OVI,emiss_CIV = region['Emission_Lya'], region['Emission_OVI'], region['Emission_CIV']

nHI,ne= region['HI_NumberDensity'], region['Electron_NumberDensity']

Z = region['Metallicity'] / 0.02


#setting the ra and dec
idx = set_central_pixel(region,pos)
dlum_0 = cosmo.luminosity_distance(zin)
dlum_0 = dlum_0.to(u.cm).value

dra = ((region['x']-pos[0])*pf['cm'])/dlum_0
ddec = ((region['y']-pos[1])*pf['cm'])/dlum_0

ra = np.zeros(len(x))+ra_0+dra
dec = np.zeros(len(y))+dec_0+ddec

#now the dlum worries about the z component
dlum = dlum_0 + z
zinterp = build_interpolator(interp_file)

z_geo = zinterp(dlum)

#Assuming v<<c
z_pec = 1. + (vlos/3.0e5)
z_app = np.zeros(len(x))+(z_cos*z_pec)


#########################################################
#### NOW I'LL WRITE OUT THE FILE ########################
#########################################################

##generate columns for binary fits file
cols = fits.ColDefs([
fits.Column(name='x', format='D', array=x, unit="Mpc/h"),
fits.Column(name='y', format='D', array=y, unit="Mpc/h"),
fits.Column(name='z', format='D', array=z, unit="Mpc/h"),
fits.Column(name='ra', format='D', array=ra, unit="radians"),
fits.Column(name='dec', format='D', array=dec, unit="radians"),
fits.Column(name='z_geo', format='D', array=z_geo, unit="1"),
fits.Column(name='z_app', format='D', array=z_app, unit="1"),
fits.Column(name='vlos', format='E', array=vlos, unit="km s^-1"),
fits.Column(name='nHI', format='E', array=nHI, unit="cm^-3"),
fits.Column(name='t', format='E', array=T, unit="K"),
fits.Column(name='metal', format='E', array=Z, unit="Zsol"),
fits.Column(name='p', format='E', array=P, unit="cgs"),
fits.Column(name='rho', format='E', array=rho, unit="g cm^-1 s^-2"),
fits.Column(name='n_el', format='D', array=ne, unit="Mpc/h"),
fits.Column(name='dlum', format='D', array=dlum, unit="Mpc/h"),
fits.Column(name='dx', format='E', array=dx, unit="Mpc/h")
fits.Column(name='emis_HI',format='D',array=emiss_HI,unit="ergs s^-1 cm^-3",
fits.Column(name='emis_OVI',format='D',array=emiss_OVI,unit="ergs s^-1 cm^-3",
fits.Column(name='emis_CIV',format='D',array=emiss_CIV,unit="ergs s^-1 cm^-3"
])



#header information for the binary fits file
ra_size = (np.max(ra)-np.min(ra))*unit['radian'].to(u.degree)
ra_size = ra_size.value
dec_size = (np.max(dec)-np.min(dec))*unit['radian'].to(u.degree)
dec_size = dec_size.value

maxdepth, mindepth = np.max(dlum),np.min(dlum) #miniumum comoving distances

## These I'm setting arbitrarily because we are using a single axis for the los
ra_los, dec_los = 0.0,0.0

## These are just in Sam's file but I don't know what they are?
ngtot, random_seed = 0.,0.

hdr = fits.Header()
hdr['RASIZE'] = ra_size
hdr['DECSIZE'] = dec_size
hdr['MAXDEPTH'] = maxdepth
hdr['MINDEPTH'] = mindepth
hdr['RANDOMSEED'] = random_seed
hdr['NGTOT'] = ngtot
hdr['RALOS'] = ra_los
hdr['DECLOS'] = dec_los



############# Actually writing the columns ####################
tbhdu = fits.BinTableHDU.from_columns(cols)
tbhdu.writeto(fileout,header=hdr,clobber=True)



