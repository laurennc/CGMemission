from lauren import *
from radial_data import *
import matplotlib.pyplot as plt

fn="/Users/laurennc/data/RyanSims/r0058_l10/redshift0058"
#fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

val = 2.37163264206e-23
pos = [0.40328598,0.47176743,0.46131516]

rad = 108.0/pf['kpc']
data = pf.h.sphere(pos,rad)

fields = ["Density","HAlphaEmissionSr","HAlphaEmissionArc","voortHalpha","HI_Number_Density"]

proj = pf.h.proj(1,fields)
width = 200./pf['kpc']
res = [1000,1000]
frb = proj.to_frb(width,res,center=pos)

xL = np.arange(-500,500)*0.2
x2 = x2[:len(x2)-1]
xL, yL = np.meshgrid(xL,xL)

rp_Arc = radial_data(frb['HAlphaEmissionArc']*(2.3529*10**(-11.)),x=xL,y=yL)
rp_Ral = radial_data(ergs_sr_TO_raleighs(frb['HAlphaEmissionSr']),x=xL,y=yL)
rp_Vor = radial_data(frb['voortHalpha'],x=xL,y=yL)

#rp_HInd = radial_data(frb['HI_Number_Density'],x=xL,y=yL)


fig,axes = plt.subplots(1,3)
thisindex = (rp_Arc.r/250. >0.01)
axes[0].plot((rp_Arc.r/250.)[thisindex],rp_Arc.mean[thisindex])
axes[0].set_yscale('log')
axes[0].set_xscale('log')
axes[0].set_xlabel('Radius / Rvir')
axes[0].set_ylabel('Surface Brightness (ergs s^-1 cm^-2 arcsec^-2)')

thisindex = (rp_Ral.r/250. >0.01)
axes[1].plot((rp_Ral.r/250.)[thisindex],rp_Ral.mean[thisindex])
axes[1].set_yscale('log')
axes[1].set_xscale('log')
axes[1].set_xlabel('Radius / Rvir')
axes[1].set_ylabel('Surface Brightness (photon s^-1 cm^-2 sr^-1)')

thisindex = (rp_Vor.r/250. >0.01)
axes[2].plot((rp_Vor.r/250.)[thisindex],rp_Vor.mean[thisindex])
axes[2].set_yscale('log')
axes[2].set_xscale('log')
axes[2].set_xlabel('Radius / Rvir')
axes[2].set_ylabel('Dense Only Surface Brightness (ergs s^-1 cm^-2 sr^-1)')


plt.savefig('radial_profile_y.png')


