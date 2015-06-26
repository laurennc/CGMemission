from yt.mods import *
import numpy as np
import cPickle
import yt.analysis_modules.halo_profiler.api as HP

fn="/u/10/l/lnc2115/vega/data/Ryan/r0054/redshift0054"
pf = load(fn, file_style="%s.grid.cpu%%04i")

#def write_halo_list(pf):
halo_list = HaloFinder(pf)
halo_list.write_out("HopAnaylsis_z02.out")
#	return

#def finding_halo_properties(fn):
#hp = HP.HaloProfiler(fn, halo_list_file='HopAnalysis.out')
#hp.add_halo_filter(HP.VirialFilter,must_be_virialized=True,
#              overdensity_field='ActualOverdensity',
#              virial_overdensity=200,
#              virial_filters=[['TotalMassMsun','>=','1e8']],
#              virial_quantities=['TotalMassMsun','RadiusMpc'])
#hp.make_profiles(filename="VirialHaloes.out")
	

