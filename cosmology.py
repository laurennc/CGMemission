import numpy as np
from astropy.cosmology import WMAP5 as cosmo
from astropy import units


def give_angular_size(z,length,unit):
	if len(length) == len(unit):
		out = [cosmo.arcsec_per_kpc_proper(z)*(length[x]*unit[x]).to(u.kpc) for x in range(len(length))]
	else:
		out = None
	return out

def give_proper_kpc(z,angle,unit):
	if len(angle) == len(unit):
		out = cosmo.kpc_proper_per_arcmin(z)*(angle[x]*unit[x]).to(u.arcmin) for x in range(len(angle))]
	else:
		out = None
	return out


