import h5py
import numpy as np

from TlacFile import TlacFile


def _join_struct_arrays(arrays):
    """Joins numpy structured arrays.

    From http://stackoverflow.com/a/5355974/1834164, much faster than
    numpy.lib.recfunctions.append_fields
    """
    sizes = np.array([a.itemsize for a in arrays])
    offsets = np.r_[0, sizes.cumsum()]
    n = len(arrays[0])
    joint = np.empty((n, offsets[-1]), dtype=np.uint8)
    for a, size, offset in zip(arrays, sizes, offsets):
        joint[:,offset:offset+size] = a.view(np.uint8).reshape(n,size)
    dtype = sum((a.dtype.descr for a in arrays), [])
    return joint.ravel().view(dtype)


class TlacHDF5(TlacFile):
    """For reading tlac hdf5 files."""

    def __init__(self, filename, dataspace):
        TlacFile.__init__(self, filename)
        self.dataspace = dataspace

        self.dat = []


    def __getitem__(self, key):
        if((self.dat == []) or (not key in self.dat.dtype.names)):
            self.load(key)
        r = self.dat[key]
        if(len(r) == 0):
            return []
        else:
            #return r[:,0]
            return r.reshape((r.shape[0],))


    def keys(self, type = 'loaded'):
        """Returns the keys (available photon data) of the instance.

        type can be `loaded` (default) to return loaded keys or
        `all` to return all available keys
        """
        if type == 'loaded':
            if(self.dat == []):
                return []
            return self.dat.dtype.names
        elif type == 'all':
            f = h5py.File(self.filename, "r")
            r = f[self.dataspace].dtype.names
            f.close()
            return r
        else:
            raise ValueError("`type` can be `loaded` or `all`.")


    def load(self, keys, force = False):
        """Not implemented yet for individual load.
        Loads all now.
        """
        if(isinstance(keys,str)):
            keys = [keys]

        # Check if already loaded
        if(self.dat != []):
            keys = [ k for k in keys if k not in self.dat.dtype.names ]
        if(keys == []):
            return
        t = tuple(keys,)

        # Load new data
        f = h5py.File(self.filename, "r")
        r = f[self.dataspace].__getitem__(t)
        if(len(t) == 1): # make to struct array if not already
            r.dtype = [ (t[0], r.dtype) ]
        f.close()

        if(self.dat == []):
            self.dat = r
        else:
            self.dat = _join_struct_arrays((self.dat, r))
        

    def load_all(self, force = False):
        if(self.dat == [] or force == True):
            f = h5py.File(self.filename, "r")
            self.dat = f[self.dataspace][:]
            f.close()

    def unload_all(self):
        self.dat = []
        

        
        
