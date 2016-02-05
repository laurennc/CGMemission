import matplotlib as mpl
mpl.use('agg')
from lauren import *
from radial_data import *
import matplotlib.pyplot as plt

#fn="/Users/laurennc/data/RyanSims/r0058_l10/redshift0058"
fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

val = 2.37163264206e-23
pos = [0.40328598,0.47176743,0.46131516]

#fields = ["HI_Number_Density","HAlpha_Emissivity","HAlpha_Emissivity_R","HAlpha_Voort_R"]
fields = ["HAlpha_Emissivity_R","HAlpha_Voort_R"]

proj = pf.h.proj(0,fields)
width = 200./pf['kpc']
res = [500,500]
frb = proj.to_frb(width,res,center=pos)

xL = np.arange(-250,250)*0.4
xL, yL = np.meshgrid(xL,xL)

#rp_Arc = radial_data(frb['HAlpha_Emissivity'],x=xL,y=yL)
rp_Ral = radial_data(frb['HAlpha_Emissivity_R'],x=xL,y=yL)
rp_Vor = radial_data(frb['HAlpha_Voort_R'],x=xL,y=yL)
#rp_hdens = radial_data(frb['HI_Number_Density'],x=xL,y=yL)


#thisindex = (rp_Arc.r/250. >0.01)
#plt.plot((rp_Arc.r/250.)[thisindex],rp_Arc.mean[thisindex])
#plt.yscale('log')
#plt.xscale('log')
#plt.xlabel('Radius / Rvir')
#plt.ylabel('Surface Brightness (ergs s^-1 cm^-2 arcsec^-2)')
#plt.savefig('rp_haArc_cloudy_smooth.png')
#plt.close()

thisindex = (rp_Ral.r/250. >0.01)
plt.plot((rp_Ral.r/250.)[thisindex],rp_Ral.mean[thisindex])
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Radius / Rvir')
plt.ylabel('Surface Brightness (R)')
plt.xlim(10**-2,1)
plt.ylim(10**-5,10)
plt.savefig('rp_haR_cloudy_smooth.png')
plt.close()

thisindex = (rp_Vor.r/250. >0.01)
plt.plot((rp_Vor.r/250.)[thisindex],rp_Vor.mean[thisindex])
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Radius / Rvir')
plt.ylabel('Surface Brightness (R)')
plt.xlim(10**-2,1)
plt.ylim(10**-5,10)
plt.savefig('rp_haVoortR_cloudy_smooth.png')
plt.close()

#thisindex = (rp_hdens.r/250. >0.01)
#plt.plot((rp_hdens.r/250.)[thisindex],rp_hdens.mean[thisindex])
#plt.yscale('log')
#plt.xscale('log')
#plt.xlabel('Radius / Rvir')
#plt.ylabel('HI Column Density (cm^-2)')
#plt.savefig('rp_hdens_cloudy_smooth.png')
#plt.close()



