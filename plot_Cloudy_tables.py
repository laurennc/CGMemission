from lauren import *
from plotting_routines import *

patt = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/bertone_exact/z000/bertone_exact_run"

#patt = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/all_lines/all_emissivity_run"
#patt =  "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/euvb/euvb_run"


#files with hden = 1, -3, -6 respectively
#inputfiles = [patt+'15.dat',patt+'7.dat',patt+'1.dat']
inputfiles = [patt+'3.dat',patt+'2.dat',patt+'1.dat']
outputfile = 'bertone_exact_double_models.png'
xlen = 3
ylen = 1

plot_Cloudy_loop(inputfiles,outputfile,xlen,ylen)


