from lauren import *
from plotting_routines import *

#patt = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/bertone_exact/bertone_exact_run"

patt = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/all_lines/all_emissivity_run"


#files with hden = 1, -3, -6 respectively
inputfiles = [patt+'15.dat',patt+'7.dat',patt+'1.dat']
outputfile = 'coronal_models.png'
xlen = 3
ylen = 1

plot_Cloudy_loop(inputfiles,outputfile,xlen,ylen)


