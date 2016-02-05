import os.path
import numpy as np

import logging

from TlacHeader import TlacHeader
from TlacHDF5 import TlacHDF5
from TlacASCII import TlacASCII

class TlacDat:
    def __init__(self, filename, verbose = False):
        """Combines tlac data, i.e., the Lya, UV data and the header (config)
        file.
        
        `filename` can be any of these files.
        """
        self.filename = os.path.expanduser(filename)
        
        self.header = []
        self.load_header()
        
        ## Set logging levels
        if(verbose):
            logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                                format='%(message)s')
            
            
        ## Some properties that can be accessed from the outside
        # Lya
        self.ndat  = float(self.header.get("result", "ndat"))
        self.nphot = float(self.header.get("run", "nphot"))
        # UV
        self.nphot_uv = float(self.header.get("run", "nphot_uv"))
        if(self.nphot_uv != 0):
            self.ndat_uv  = float(self.header.get("result", "ndat_uv"))
        else:
            self.ndat_uv  = 0

        ## Bugfix: An older version of tlac writes nphot < 0 although
        ##         zero photons are written
        self.ndat     = np.max([self.ndat, 0])
        self.nphot    = np.max([self.nphot, 0])
        self.ndat_uv  = np.max([self.ndat_uv, 0])
        self.nphot_uv = np.max([self.nphot_uv, 0])
            
        ## Connect data
        self.dat    = []
        self.dat_uv = []
        self.dat = self.set_dat(self.header.get("output", "lya"))
        if(self.ndat_uv > 0):
            self.dat_uv = self.set_dat(self.header.get("output", "uv"))

            
    def __getitem__(self, key):
        """Access a property "X" of "Lya" either via TlacDat["Lya","X"] or
        TlacDat["Lya"]["X"]
        """
        if(isinstance(key, str)):
            if(key == 'Lya'):
                return self.dat
            elif(key == 'UV'):
                return self.dat_uv
            else:
                raise Exception("Unrecognized key %s") % key
        elif(len(key) == 2):
            return self.__getitem__(key[0])[key[1]]
        else:
            raise Exception("Too many keys!")
                

    def keys(self, **kwargs):
        """Returns keys.
        """
        return self.dat.keys(**kwargs)
            
    def set_dat(self, key):
        """Sets dat type according to config file.
        """
        filepre = ".".join(self.filename.split(".")[:-1])
        
        key = key.upper()
        
        if(key == "ASCII_LYA"):
            to_set = TlacASCII( filepre + ".dat" )
        elif(key == "ASCII_UV"):
            to_set = TlacASCII( filepre + ".uvdat" )
        elif(key == "HDF5_LYA"):
            to_set = TlacHDF5( filepre + ".h5", "Lya")
        elif(key == "HDF5_UV"):
            to_set = TlacHDF5( filepre + ".h5", "UV" )
        else:
            raise Exception("Key "+ key + " not recognized!")
            
        return to_set
            
    def load_header(self):
        """Loads file header.
        """
        tmp = self.filename.split(".")
        cfgfn = ".".join(tmp[:-1]) + ".cfg"
        
        self.header = TlacHeader(cfgfn)
        
        
    def compare_with(self, to_compare):
        """Compares the header (parameters) of the data set with
        the one of `to_compare`.
        
        Returns a list of the keys that are not equal.
        """
        
        not_equal = []
        
        for k, v in self.header.iteritems():
            if(to_compare.header[k] != v):
                not_equal.append(k)
                
                
        return not_equal

        
    def get_nlines(self, comment_pre = "#"):
        """Counts lines that dont start with `comment_pre` and returns result.
        """
        f = open(self.filename, "r+")
        buf = mmap.mmap(f.fileno(), 0)
        lines = 0
        while True:
            l = buf.readline()
            if(l == ''):
                break
            if(l[0:len(comment_pre)] != comment_pre):
                lines += 1
        return lines


    def load(self, phot_type = 'both', *args):
        if(phot_type == 'both'):
            self.dat.load(*wargs)
            self.dat_uv.load(*args)
        elif(phot_type == 'Lya'):
            self.dat.load(*args)
        elif(phot_type == 'UV'):
            self.dat_uv.load(*args)
        else:
            raise TypeError("Choose 'both', 'Lya' or 'UV' for phot_type.")


    def load_all(self, phot_type = 'both'):
        """Loads all photons of `phot_type`.
        phot_type can be 'both', 'Lya' or 'UV'.
        """
        if(phot_type == 'both'):
            self.dat.load_all()
            self.dat_uv.load_all()
        elif(phot_type == 'Lya'):
            self.dat.load_all()
        elif(phot_type == 'UV'):
            self.dat_uv.load_all()
        else:
            raise TypeError("Choose 'both', 'Lya' or 'UV' for phot_type.")
            
            
    def unload_all(self):
        if(self.dat != []):
            self.dat.unload_all()
        if(self.dat_uv != []):
            self.dat_uv.unload_all()


        
