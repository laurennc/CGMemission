import pymses
import numpy as np
import math
import matplotlib.pyplot as plt
from pymses.sources.ramses.output import *
from pymses.filters import RegionFilter
from pymses.filters import CellsToPoints
from pymses.filters import PointFunctionFilter
from pymses.analysis import sample_points
from pymses.utils.regions import Sphere
from pymses.analysis.visualization import *
from pymses.utils import constants as C
from astropy.io import fits
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.units as u
from astropy.io import fits
from scipy import interpolate

from classTable import *
from matplotlib.mlab import griddata
import cPickle

from readHalo2 import *

filename_cat = 'catalog.fits'

kBoltzmann = const.k_B.cgs.value ### 1.3806488*math.pow(10,-16) ### erg/K, CGS units
hydrogenMass = 1.6605402*math.pow(10,-24) ### g, CGS units

### from Grevesse+10
Xsol = 0.7380
Ysol = 0.2484
Zsol = 0.0134

############
############
############
############

### los => choose one axis to look at, e.g. x = (1,0,0)
los = (1,0,0)
RA = (0,1,0)
DEC = (0,0,1)

### center of the box at redshift defined by redshift of the output
z_center = 1.22 

############
############
############
############

# Read additional tracers : metallicity
RamsesOutput.amr_field_descrs_by_file = \
{   "3D": {"hydro" : [ Scalar("rho", 0), Vector("vel", [1, 2, 3]), Scalar("P", 4), Scalar("Z", 5) ],\
         "grav"  : [ Vector("g", [0, 1, 2]) ]
        }
}

print  '\n------------------------------\n    Loading Ramses output\n--------------------------------\n'
output = 23
ro = pymses.RamsesOutput("/galex_bingo/Simulations/boxlen100_n128_lcdmw5_nonthermal_lessref_zoom-bingo-12/", output)
info = ro.info

cosmo = FlatLambdaCDM(H0=info["H0"], Om0=info["omega_m"], Ob0=info["omega_b"])

############
############
############
############


pathToHalo = '/galex_bingo/Simulations/boxlen100_n128_lcdmw5_nonthermal_lessref_zoom-bingo-12/Halos/23/haloProps.023'
haloList = np.genfromtxt(pathToHalo, skip_header = 1)

halos = []

for i in range(len(haloList[:,0])):
	h1 = haloList[i,:]	
	h = Halo(h1[0],h1[1],h1[2],h1[3],h1[4],h1[5],h1[6],h1[7],h1[8],h1[9],h1[10],h1[11],h1[12],h1[13],h1[14])
	halos.append(h)

print '%i halos created.'%len(halos)

maxMass = -1
for h in halos:
        if maxMass < h.mass:
                maxMass = h.mass
                massiveHalo = h

halo = massiveHalo

print '\nCreating sphere region from more massive Halo'
center = np.array([halo.x/info["unit_length"].express(C.Mpc)+.5,halo.y/info["unit_length"].express(C.Mpc)+.5,halo.z/info["unit_length"].express(C.Mpc)+.5])
print 'center: %f ,%f, %f'%(center[0],center[1],center[2])
radius = halo.virialRadius(z_center)/info["unit_length"].express(C.cm)
print 'Virial radius: %f kpc    or     %f in box units'%(radius * info["unit_length"].express(C.kpc), radius)
region = Sphere(center, 3*radius)



############
############
############
############



log = np.vectorize(math.log10)
sqrt = np.vectorize(np.sqrt)
atan = np.vectorize(math.atan)
asin = np.vectorize(math.asin)

def derive_Pressure(dset):
       return  log(dset["P"]*info["unit_pressure"].express(C.kg/C.m/C.s/C.s))

def derive_Temperature(dset):
        P_cgs = dset["P"]*info["unit_pressure"].express(C.g/C.cm/C.s/C.s) ## cgs units
        rho_cgs = dset["rho"]* info["unit_density"].express(C.g/C.cm**3) ### cgs_units

        ToverMu = (hydrogenMass/kBoltzmann)*P_cgs/rho_cgs

        return log(ToverMu)

def derive_nH(dset):
        rho_cgs = dset["rho"]* info["unit_density"].express(C.g/C.cm**3) ### cgs_units
        nH = (Xsol/hydrogenMass)*rho_cgs

        return log(nH)

def derive_rho(dset):
        rho_cgs = dset["rho"]* info["unit_density"].express(C.g/C.cm**3) ### cgs_units

        return log(rho_cgs)

def derive_metallicity(dset):
        Z = dset["Z"]
        return log(Z/Zsol)



############
############
############
############


print '\nSelecting amr fields'
amr = ro.amr_source(["rho", "vel", "P", "Z"])





# AMR data filtering
print '\nFiltering the data with region'
amr_region = RegionFilter(region, amr)

#####
cell_source = CellsToPoints(amr_region)
cells = cell_source.flatten()

print 'number of cells in selected region: %i'%cells.npoints
#print cells.get_sizes()


############
############
############
############





### log density
print 'deriving log rho...'
rho_cell = derive_rho(cells)

### log HI density
print 'deriving log nH...'
HI_density_cell = derive_nH(cells)

### log temperature
print 'deriving log T/mu...'
temp_cell = derive_Temperature(cells)

### log metallicity
print 'deriving log Z/Zsol...'
metallicity_cell = derive_metallicity(cells)

### log pressure
print 'deriving log P'
pressure_cell = derive_Pressure(cells)





############ verifier que les X_, Y_, Z_ sont bien en Mpc/h!!!!!


### X_cone, Y_cone and Z_cone
print 'deriving points coordinates...'
cells_points_Mpc = (cells.points-center) * info["unit_length"].express(C.Mpc)

### conversion boite proper ou comoving?



### luminosity distance
print 'deriving comoving distances...'

## The redshift corresponds to the center point of the simulation box
X_cone_center = cosmo.comoving_distance(z_center).value ## Mpc



#AMR_position = (xcell,ycell,zcell)
print 'deriving AMR coordinates...'

### passer en matriciel
X_cone = np.array([np.dot(los,a) for a in cells_points_Mpc]) + np.array(X_cone_center)
Y_cone = np.array([np.dot(RA,a) for a in cells_points_Mpc])
Z_cone =np.array([np.dot(DEC,a) for a in cells_points_Mpc])







### zgeo

#### tabuler les z
## appeler 2 fois z_at_value pour avoir le zmin et zmax
## tabuler les distance en fonction de z, assez fin
## faire du linearInterp


d_cone = sqrt(X_cone**2 + Y_cone**2 + Z_cone**2) ## comoving distance to each cell, Mpc

def comoving2zgeo(d):
        return 1.

zgeo_cell = X_cone/X_cone #comoving2zgeo(d_cone)

dL_cell = d_cone * (1.+zgeo_cell)





### zapp 
print 'deriving zapp...'

vcell_boxunits = cells["vel"]
v_los_boxunits = np.array([np.dot(a,los) for a in vcell_boxunits])
v_los_kms = v_los_boxunits * (info["unit_length"].express(C.km)/info["unit_time"].express(C.s))

c_kms = C.c.express(C.km/C.s)

#zdopp_cell = np.array([math.sqrt((1.+v/c_kms)/(1.-v/c_kms))-1 for v in v_los_kms])
zdopp_cell = np.array([1.+v/c_kms for v in v_los_kms])

zapp_cell = zgeo_cell * zdopp_cell





### Ra, Dec

print 'deriving Ra, dec...'
#arcsec_per_kpc_proper = np.array([cosmo.arcsec_per_kpc_proper(a).value for a in zgeo_cell]) ### arcsec / kpc    ###############################################
#radian_per_Mpc_proper = arcsec_per_kpc_proper * u.arcsec.to(u.radian) * u.kpc.to(u.Mpc) 

#radian_per_Mpc_proper = cosmo.arcsec_per_kpc_proper(z_center).value * u.arcsec.to(u.radian) / u.kpc.to(u.Mpc) 

## in flat Universe, angles are the same

ra_cell = atan(Y_cone/X_cone) #Y_cone * radian_per_Mpc_proper
dec_cell = asin(Z_cone/X_cone) #Z_cone * radian_per_Mpc_proper






### cell size
print 'deriving cell sizes...'
size_cell = cells.get_sizes() * info["unit_length"].express(C.Mpc)


















###############
###############
###############    rajouter la colonne n_el => a tabuler par Cloudy
###############
###############

n_el_cell = 0 * X_cone


### generate binary fits file


cols = fits.ColDefs([
        fits.Column(name='x', format='D', array=X_cone, unit="Mpc/h"),
fits.Column(name='y', format='D', array=Y_cone, unit="Mpc/h"),
fits.Column(name='z', format='D', array=Z_cone, unit="Mpc/h"),
fits.Column(name='ra', format='D', array=ra_cell, unit="radians"),
fits.Column(name='dec', format='D', array=dec_cell, unit="radians"),
fits.Column(name='z_geo', format='D', array=zgeo_cell, unit="1"),
fits.Column(name='z_app', format='D', array=zapp_cell, unit="1"),
fits.Column(name='vlos', format='E', array=v_los_kms, unit="km s^-1"),
fits.Column(name='nHI', format='E', array=HI_density_cell, unit="cm^-3"),
fits.Column(name='t', format='E', array=temp_cell, unit="K"),
fits.Column(name='metal', format='E', array=metallicity_cell, unit="Zsol"),
fits.Column(name='p', format='E', array=pressure_cell, unit="cgs"),
fits.Column(name='rho', format='E', array=rho_cell, unit="g cm^-1 s^-2"),
fits.Column(name='n_el', format='D', array=n_el_cell, unit="Mpc/h"),
fits.Column(name='dlum', format='D', array=dL_cell, unit="Mpc/h"),
fits.Column(name='dx', format='E', array=size_cell, unit="Mpc/h")
])




############
############
############
############

#emissivities

tablePath = 'tables/'
list_lines = [line.strip() for line in open(tablePath + 'listLines.lines')]

idLine = 0 ## TODO

print 'Reading emissivity table %s...'%list_lines[idLine]
tab = [line.strip() for line in open(tablePath + 'emissivityTable_%s.tab'%list_lines[idLine])]


n_info=tab[0].split(' ')
n_start = float(n_info[0])
n_step = float(n_info[1])
n_nb = float(n_info[2])

n = np.array(np.linspace(n_start,n_start+(n_nb-1)*n_step,n_nb))

t_info=tab[1].split(' ')
t_start = float(t_info[0])
t_step = float(t_info[1])
t_nb = float(t_info[2])

t = np.array(np.linspace(t_start,t_start+(t_nb-1)*t_step,t_nb))


em = []
for i in range(len(n)):
    em.append([float(a) for a in tab[i+2].split(' ')])


points = []
values = []
for i in range(len(n)):
        for j in range(len(t)):
                points.append([n[i],t[j]])
                values.append(em[i][j])



print 'interpolating %s emissivity on cells...'%list_lines[idLine]
f = interpolate.LinearNDInterpolator(points, values, fill_value = -60)


print 'deriving emissivity %s...'%list_lines[idLine]
emissivity = f([[a,b] for (a,b) in zip(HI_density_cell, temp_cell)])


### add emissivities to catalog
cols.add_col(fits.Column(name='emissivity_%s'%list_lines[idLine], format='D', array=emissivity, unit="erg cm^-3 s^-1"))










plt.figure()
plt.hist(emissivity, bins=500, histtype='step', normed=1)
plt.title('emissivity distribution %s'%list_lines[idLine])
plt.xlabel('em(%s)'%list_lines[idLine])

plt.figure()
plt.plot(temp_cell, HI_density_cell,'+')
plt.xlim([3.,8.])
plt.ylim([-7.,2.])
plt.ylabel('log nHI')
plt.xlabel('log T/mu')
plt.title('phase diagram')


### emissivity table


plt.imshow(em, interpolation = 'nearest', vmin = -45., extent=[3.,8.,-7.,2.])
plt.colorbar()
#plt.title("log em(%s)"%table._lineNames[idLine])

#plt.xlabel('log T')
#plt.ylabel('log nH')

#plt.figure()
#plt.plot(temp_cell, HI_density_cell,'+')

### metallicity distribution

plt.figure()
plt.hist(metallicity_cell,bins=500)
plt.xlabel('log Z/Zsol')
plt.title('metallicity distribution')
plt.show()

############
############
############
############





tbhdu = fits.BinTableHDU.from_columns(cols)
tbhdu.writeto(filename_cat, clobber=True)

print 'catalog created.'




###########
##########
###########

### Read catalog
print 'Reading catalog...'

hdulist = fits.open(filename_cat)
tbdata = hdulist[1].data
cols = hdulist[1].columns
print cols.info()


### add keywords

ra_size = 0.
dec_size = 0.
maxdepth = 0.
mindepth = 0.
random_seed = 0.
ngtot = 0.
ra_los = 0.
dec_los = 0.

hdulist[1].header['ra_size'] = ra_size
hdulist[1].header['dec_size'] = dec_size
hdulist[1].header['maxdepth'] = maxdepth
hdulist[1].header['mindepth'] = mindepth
hdulist[1].header['random_seed'] = random_seed
hdulist[1].header['ngtot'] = ngtot
hdulist[1].header['ra_los'] = ra_los
hdulist[1].header['dec_los'] = dec_los



fits.writeto(filename_cat, hdulist[1].data, hdulist[1].header, clobber=True)
