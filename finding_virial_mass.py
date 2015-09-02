import matplotlib
matplotlib.use('Agg')
from yt.mods import *
import numpy as np
import cPickle
import yt.analysis_modules.halo_profiler.api as HP
from yt.analysis_modules.halo_finding.api import *

##################### z = 0 #########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0058_l10/redshift0058"
#pos = [0.40328598,0.47176743,0.46131516]
##################### z = 0.2 ##########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
####################  z = 0.5 ###########################
#fn="/u/10/l/lnc2115/vega/data/Ryan/r0048/redshift0048"
####################  z = 1.0 ##########################
fn="/u/10/l/lnc2115/vega/data/Ryan/r0038/redshift0038"

pf = load(fn, file_style="%s.grid.cpu%%04i")

#def write_halo_list(pf):
halo_list = HaloFinder(pf)
#halo_list.write_out("HopAnaylsis_z1.out")
halo_list.dump(basename='HopAnalysis_z1')
#	return

#halos = LoadHaloes(pf,'HopAnalysis_z0')

#def finding_halo_properties(fn):
#hp = HP.HaloProfiler(fn, halo_list_file='HopAnalysis.out')
#hp.add_halo_filter(HP.VirialFilter,must_be_virialized=True,
#              overdensity_field='ActualOverdensity',
#              virial_overdensity=200,
#              virial_filters=[['TotalMassMsun','>=','1e8']],
#              virial_quantities=['TotalMassMsun','RadiusMpc'])
#hp.make_profiles(filename="VirialHaloes.out")
	

