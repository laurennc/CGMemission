import matplotlib as mpl
mpl.use('agg')
from lauren import *
from radial_data import *
import matplotlib.pyplot as plt

fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

val = 2.37163264206e-23
pos = [0.40328598,0.47176743,0.46131516]

#fields = ["HI_Number_Density","HAlpha_Emissivity","HAlpha_Emissivity_R","HAlpha_Voort_R"]
fields = ["HAlpha_Emissivity_R","HAlpha_Voort_R"]
width = 200./pf['kpc']
#res = [500,500]
#res =[ 200,200]
res = [40,40]

xL = np.arange(-20,20)*10.0
#xL = np.arange(-100,100)
#xL = np.arange(-250,250)*0.4
xL, yL = np.meshgrid(xL,xL)
r = abs(xL+1j*yL)

projx = pf.h.proj(0,fields)
frbx = projx.to_frb(width,res,center=pos)
rp_Ralx = radial_data(frbx['HAlpha_Emissivity_R'],x=xL,y=yL)
rp_Vorx = radial_data(frbx['HAlpha_Voort_R'],x=xL,y=yL)

projy = pf.h.proj(1,fields)
frby = projy.to_frb(width,res,center=pos)
rp_Raly = radial_data(frby['HAlpha_Emissivity_R'],x=xL,y=yL)
rp_Vory = radial_data(frby['HAlpha_Voort_R'],x=xL,y=yL)

projz = pf.h.proj(2,fields)
frbz = projz.to_frb(width,res,center=pos)
rp_Ralz = radial_data(frbz['HAlpha_Emissivity_R'],x=xL,y=yL)
rp_Vorz = radial_data(frbz['HAlpha_Voort_R'],x=xL,y=yL)


rp_Ralmean = (rp_Ralx.mean + rp_Raly.mean + rp_Ralz.mean)/3.0
rp_Vormean = (rp_Vorx.mean + rp_Vory.mean + rp_Vorz.mean)/3.0


#thisindex = (rp_Ralx.r/250. >0.01)
#plt.plot((rp_Ralx.r/250.)[thisindex],rp_Ralmean[thisindex])
plt.plot(rp_Ralx.r,rp_Ralmean)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Radius (kpc)')
plt.ylabel('Surface Brightness (R)')
#plt.xlim(10**-2,1)
plt.xlim(0,200)
plt.ylim(10**-5,10)
plt.savefig('rp_haR_cloudy_avg_smooth10.png')
plt.close()

#thisindex = (rp_Vorx.r/250. >0.01)
#plt.plot((rp_Vorx.r/250.)[thisindex],rp_Vormean[thisindex])
plt.plot(rp_Vorx.r,rp_Vormean)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Radius (kpc) ')
plt.ylabel('Surface Brightness (R)')
plt.xlim(0,200)
#plt.xlim(10**-2,1)
plt.ylim(10**-5,10)
plt.savefig('rp_haVoortR_cloudy_avg_smooth10.png')
plt.close()




