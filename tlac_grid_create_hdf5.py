"""Functions for creating a tlac_grid in hdf5 format.
"""

import time, sys, os
import numpy as np
import h5py

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                 '../../')))
import tlac_analysis as ta    

from tlac_grid_create import *


def open_file(fn):
    """Opens hdf5 file `fn` and returns pointer
    """
    return h5py.File(fn, "w")
    

def write_header(f, comment = None):
    """Writes full header (with blocks) to file `f`. If no comment is given
    takes current time.
    """
    if comment is None:
        comment = "tlac-grid created with python script on %s\0" % time.ctime()

    f.create_dataset("header", (0,))
    dset = f['header']
    
    dset.attrs["comment"] = np.string_(comment)
    dset.attrs["tlac_o_grid_magic_number"] = np.int32(tlac_magicnumber())
    
    
def write_grid_data(f, grid_type, ncells, boxsize, L_tot):
    """Writes full grid data (with blocks) to file `f`.

    boxsize should be given in cm
    """
    f.create_dataset("grid_data", (0,))
    dset = f['grid_data']
    
    dset.attrs["grid_type"] = np.int32(grid_type)
    dset.attrs["ncells"] = [ np.int32(i) for i in ncells ]
    dset.attrs["boxsize"] = [ float(i) for i in boxsize ]
    dset.attrs["L_tot"] = float(L_tot)
    

def write_cell_data(f, pos_coords, L, T, n_HI, n_D, vel, n_H, Z, nclouds = None,
                    verbose = True, dset_dict = {'compression' : 4 }):
    """Writes full cell data section to file `f`.

    Keyword arguments:
    f                    --  File handler to write in (hdf5)
    L, T, n_HI, n_D, vel --  numpy.arrays with shape (N1, N2, N3)
                             ((N1,N2,N3,3) for `vel`).
                             `T` in K, others in cgs (`vel` will be converted)
                             Note that tlac hdf5 grid supports not storing a
                             certain value (must be given) by tlac directly
                             then. For not storing something, set it to `None`.
    nclouds               -- Number of clouds in the cells.
    pos_coords            -- List of 1D arrays with lengths N1, N2, N3. Center
                             of cell along this axis.
    dset_dict             -- Will be passed to `hyp5.create_dataset`. Can, e.g.,
                             contain compression filters. Default is moderate
                             compression and fast speed.
    """
    shape = np.array([len(i) for i in pos_coords])
    ntot = np.product(shape)

    ### fix here if slow!
    pos_grid = np.zeros(np.concatenate( (shape, [3]) ))
    pgrid = np.meshgrid(pos_coords[0], pos_coords[1], pos_coords[2],
                        indexing = 'ij')
    pos_grid[:,:,:,0], pos_grid[:,:,:,1], pos_grid[:,:,:,2] = pgrid
    f.create_dataset("cell_data/pos", data = pos_grid,  **dset_dict)
    
    if L is not None:
        f.create_dataset("cell_data/L", data = L, **dset_dict)
    if T is not None:
        f.create_dataset("cell_data/T", data = T, **dset_dict)
    if n_HI is not None:
        f.create_dataset("cell_data/n_HI", data = n_HI, **dset_dict)
    if n_D is not None:
        f.create_dataset("cell_data/n_D", data = n_D, **dset_dict) 
    if vel is not None:
        v_th = ta.thermal_velocity(T)
        vel[:,:,:,0] /= v_th
        vel[:,:,:,1] /= v_th
        vel[:,:,:,2] /= v_th
        f.create_dataset("cell_data/velocity", data = vel, **dset_dict)
    if nclouds is not None:
        f.create_dataset("cell_data/nclouds", data = nclouds.astype(int),
                         **dset_dict)
    ####LAUREN ADDED HERE####
    if n_H is not None:
        f.create_dataset("cell_data/n_H",data=n_H, **dset_dict)
    
    if Z is not None:
        f.create_dataset("cell_data/Z",data=Z, **dset_dict)


def write_cloud_data(f):
    """Does nothing except empty block so far...    
    """
    
def close_file(f):
    f.close()
 
