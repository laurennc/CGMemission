from yt.mods import *
import numpy as np
import cPickle

fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
pf = load(fn, file_style="%s.grid.cpu%%04i")


halo_list = HaloFinder(pf)
halo_list.write_out("HopAnaylsis.out")


