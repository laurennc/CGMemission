import matplotlib
matplotlib.use('Agg')
from lauren import *
import matplotlib.pyplot as plt

fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

pos = [0.40328598,0.47176743,0.46131516]
rad = 108.0/pf['kpc']

data_source = pf.h.sphere(pos,rad)

plt.scatter(np.log10(data_source['Density']),np.log10(data_source['Temperature']),c=np.log10(data_source['CellMassMsun']),marker='.',lw=0,alpha=0.25)

plt.xlabel('log(Density)')
plt.ylabel('log(Temperature)')
cb = plt.colorbar()
cb.set_label('log(Cell Mass) (M_sun)')
plt.savefig('phase_plot.png')
plt.close()



