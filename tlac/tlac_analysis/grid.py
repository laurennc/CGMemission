"""File for functions analyzing a tlac grid file.
"""

from physics import thermal_velocity
import numpy as np
import h5py

from prettytable import PrettyTable


def grid_load(filename):
    """
    Reads tlac grid file.
    
    Keyword Arguments:
    filename -- Name of grid file to load

    Returns:
    Dictionary of the following form:
    {
      'cell_dat' : {'T', 'vel', ...} (which are 1D and 3D arrays)
    }

    Notes:
    This function will call either the ascii or the hdf5 read version.
    The common format described above is not realized yet. Depending on the
    file format the resulting data is stored in 1D (ascii) or 3D (hdf5) arrays.
    If this file is needed more often, one should maybe find a common ground...
    """
    ext = filename.split(".")[-1]
    if ext == "grid" or ext == "dat":
        return _grid_load_ascii(filename)
    elif ext == "hdf5" or ext == "h5":
        return _grid_load_hdf5(filename)

    print "[Warning] Could not recognize grid filetype based on filename. "\
        "Trying now ascii..."

    try:
        return _grid_load_ascii(filename)
    except:
        print "[Warning] Loading as ascii failed! Trying now hdf5"

    try:
        return _grid_load_hdf5(filename)
    except:
        print "[Warning] Loading as hdf5 failed!"

    raise Exception("Grid \"%s\" could not be loaded." %(filename))        


def _grid_load_ascii(filename):
    f = open(filename,"r")
    lines = f.readlines()
    f.close()

    cell_start = lines.index("### START CELL_DATA ###\n")
    cell_end = lines.index("### END CELL_DATA ###\n")

    celldat = np.array([np.fromstring(i, sep=' ')
                        for i in lines[(cell_start+1):cell_end] ])

    cell_dict = {}
    cell_dict['L'] = celldat[:,0]
    cell_dict['T'] = celldat[:,1]
    cell_dict['n_HI'] = celldat[:,2]
    cell_dict['n_D'] = celldat[:,3]
    cell_dict['pos'] = celldat[:,4:7]
    cell_dict['velocity'] = celldat[:,7:10]
    cell_dict['nclouds'] = celldat[:,10]

    # some extra info
    cell_dict['vel_conv'] = np.array([cell_dict['velocity'][:,i] *
                                      thermal_velocity(cell_dict['T'])
                                      for i in range(3) ]) 
    

    return {
        'cell_dat' : cell_dict
    }


def _grid_load_hdf5(fn):
    f = h5py.File(fn, "r")

    cell_dict = {}
    for name, dset in f['cell_data/'].iteritems():
        cell_dict[name] = dset[...]
    f.close()

    if 'velocity' in cell_dict and 'T' in cell_dict:
        cell_dict['vel_conv'] = np.array([cell_dict['velocity'][:,:,:,i] *
                                          thermal_velocity(cell_dict['T'])
                                          for i in range(3) ]) 

    return {
        'cell_dat' : cell_dict
    }


def velocity_rms(velocity):
    """ Returns 3d rms of velocity data in input units
    
    Argument:
    velocity -- Array of velocity data
    """

    return np.sqrt(3 * np.mean(velocity**2))


def grid_summary(griddat):
    """
    Prints tlac grid summary.
    
    Keyword Arguments:
    griddat -- Grid data in dictionary form as returned from
               `tlac_analysis.grid_load`
    """
    print "Key\tMean\tstd\tmin\tmax"
    print "------------------------------------------------------------"
    for k, v in griddat['cell_dat'].iteritems():
        print "%s\t%.3e\t%.3e\t%.3e\t%.3e" %(k, np.mean(v), np.std(v),
                                             np.min(v), np.max(v))
    print "------------------------------------------------------------"
    print "Some other things:"
    print "\tv_rms (km/s)",velocity_rms(griddat['cell_dat']['vel_conv']) / 1e5

    
    

    
