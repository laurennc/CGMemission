import matplotlib as mpl
mpl.use('agg')
from lauren import *
from radial_data_lauren import *
import matplotlib.pyplot as plt

fn="/hpc/astrostats/astro/users/lnc2115/Ryan/r0058_l10/redshift0058"
pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

val = 2.37163264206e-23
pos = [0.40328598,0.47176743,0.46131516]

#fields = ["HI_Number_Density","HAlpha_Emissivity","HAlpha_Emissivity_R","HAlpha_Voort_R"]
fields = ["HAlpha_Emissivity_R","HAlpha_Voort_R"]
width = 200./pf['kpc']
#res = [500,500]
res = [200,200]

xL = np.arange(-100,100)
xL, yL = np.meshgrid(xL,xL)
r = abs(xL+1j*yL)

proj = pf.h.proj(0,fields)
frb = proj.to_frb(width,res,center=pos)
rp_Ralx = radial_data(frb['HAlpha_Emissivity_R'],x=xL,y=yL)
rp_Vorx = radial_data(frb['HAlpha_Voort_R'],x=xL,y=yL)

proj = pf.h.proj(1,fields)
frb = proj.to_frb(width,res,center=pos)
rp_Raly = radial_data(frb['HAlpha_Emissivity_R'],x=xL,y=yL)
rp_Vory = radial_data(frb['HAlpha_Voort_R'],x=xL,y=yL)

proj = pf.h.proj(2,fields)
frb = proj.to_frb(width,res,center=pos)
rp_Ralz = radial_data(frb['HAlpha_Emissivity_R'],x=xL,y=yL)
rp_Vorz = radial_data(frb['HAlpha_Voort_R'],x=xL,y=yL)


rp_Ralmean = (rp_Ralx.mean + rp_Raly.mean + rp_Ralz.mean)/3.0
rp_Vormean = (rp_Vorx.mean + rp_Vory.mean + rp_Vorz.mean)/3.0


#CAN ALSO OVERPLOT THE MEAN TREND FROM ABOVE
#LINESTYLE='-', plt.set_linewidth(1.2)

rflat = r.flatten()
rHalpha = frb['HAlpha_Emissivity_R'].flatten()
#rpRalxflat = rp_Ralx.flatten()
thisindex = (rflat/250. > 0.01)
thisindex2 = (rp_Ralx.r/250. > 0.01)
print thisindex2
plt.plot(rflat,rHalpha,'bo',markersize=0.30)
plt.plot(rp_Ralx.r,rp_Ralmean,linewidth=1.2)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Radius (kpc)')
plt.ylabel('Surface Brightness (R)')
plt.axis([1,200,10**-5,10])
#plt.xlim(10**0.,10**2.3)
#plt.xlim(10**-2,1)
#plt.ylim(10**-5,10)
plt.savefig('rp_radialscatter_x.png')
plt.close()
