import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np

ions = ['Lya','C3_977','CIV_1548','OVI_1032']
colnum = [1,7,3,5]
hden = [-5,-3,-2]

pattbeg = '/u/10/l/lnc2115/vega/data/Ryan/cloudy_out/grid_galquas/emis/z0/g1q1'
pattmid = '/g1q1_run9.dat'
inputfile = pattbeg+pattmid

data = np.genfromtxt(inputfile,skip_header=12)
x = data[:,0]

plt.plot(x,data[:,1],linewidth=1.5)

plt.xlim(3,7)
plt.savefig('sam_compare.png')
plt.close()

