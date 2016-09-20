import matplotlib
matplotlib.use('Agg')
from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt

fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
pf = load(fn, file_style="%s.grid.cpu%%04i")

pos = [0.40320258,0.4718025,0.46123817]
rad = 320.

sp = pf.h.sphere(pos,rad/pf['kpc'])

pc = PlotCollection(pf,pos)
nhT = pc.add_phase_sphere(rad,'kpc',['H_NumberDensity','Temperature','cooling_time'],weight=None,x_bins=256,y_bins=256)

nhT.set_xlim(10**-6.,10**2.5)
nhT.set_ylim(10**3.75,10**6.8)
nhT.set_zlim(10**3.,10**20.)

pc.save('phaseplot_coolingtime_smallcmap')




