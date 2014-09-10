#from lauren import *
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from plotting_routines import *

#patt = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/bertone_exact/z000/bertone_exact_run"

#pattbeg = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/Ions/"
#pattmid = ["euvbIon_factor001/","euvbIon_factor1/","control/bertone/","euvbIon_factor10/","euvbIon_factor1000/"]
#pattend = "euvb_ion_run"
#pattmid = ["bertIon_factor001/","bertIon_factor1/","bertIon_control/","bertIon_factor10/","bertIon_factor1000/"]
pattend = "bert_ion_run"
pattmid = ["holder/","HMgals/","HMquas/"]
pattbeg = "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/Ions/varyShape/"

ion = 'C'
outputfile = 'HM05_ionfrac_'+ion+'.png'
#outputfile = 'bert_ionfrac_'+ion+'.png'

#files with hden = 1, -3, -6 respectively
#inputfiles = [patt+'15_'+ion+'.dat',patt+'7_'+ion+'.dat',patt+'1_'+ion+'.dat']
#inputfiles = [patt+'3.dat',patt+'2.dat',patt+'1.dat']
inputPatts = [pattbeg+x+pattend for x in pattmid]
inputPatts[0] =  "/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/Ions/bertIon_control/bert_ion_run"
#inputPatts[2] = pattbeg+pattmid[2]+"bert_ion_run"
xlen = 3
ylen = 1
#ylen = 5

#plot_Cloudy_loop(inputfiles,outputfile,xlen,ylen)
#outputfile = "varyShape_"+ion+".png"
plot_ion_loop(inputPatts,outputfile,ion,xlen,ylen,special=True)

