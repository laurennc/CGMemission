
from yt.mods import *

fn="r0058_l10/redshift0058"

pf = load(fn, file_style="%s.grid.cpu%%04i") # load data

SlicePlot(pf, 'x', "Density", width = (800.0, 'kpc')).save()
SlicePlot(pf, 'y', "Density", width = (800.0, 'kpc')).save()
SlicePlot(pf, 'z', "Density", width = (800.0, 'kpc')).save()
