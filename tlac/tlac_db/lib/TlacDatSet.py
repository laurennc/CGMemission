## @class TlacDatSet
## A set of tlac data files containing all the same tlac run parameters.
##
import fileinput
import glob
import sys
import numpy as np

import os.path
import os, sys, inspect

from TlacDat import TlacDat

class TlacDatSet:
    def __init__(self, file_prefix = "", filelst = [], verbose = False):
        """ Inits tlac set and adds files starting with `file_prefix`
        (can be, e.g., directory).

        Alternatively, `filelst` can be given to add files directly.
        
        `verbose` can be set to true if printout is wanted
        """
        # Individual file objects
        self.tlac_dat = []
        
        ## Header object
        self.header = []

        ## Print out stuff
        self.verbose = verbose

        if(file_prefix):
            self.add_files(glob.glob(file_prefix + "*.cfg"))

        if(filelst):
            self.add_files(filelst)


    def __getitem__(self, (phot_type, prop)):
        return np.concatenate([ i[phot_type][prop] for i in self.tlac_dat ])

    
    def keys(self, **kwargs):
        """Returns list of photon properties that are to obtain
        """
        return self.tlac_dat[-1].keys(**kwargs)
        
            
    def add_files(self, filelst, **kwargs):
        """Adds several tlac  files.
        Checks if parameter settings are the same and loads UV data if existing.
        """
        if(isinstance(filelst, basestring)):
            filelst = [filelst]
        
        [ self.add_file(i, **kwargs) for i in filelst ]


    def add_file(self, filename, check = True):
        """Adds file to TlacDatSet.
        Returns True if added otherwise False

        Checks if ok to add (same parameters) iff `check=True` (default)
        """
        tdl = self.tlac_dat

        if(self.verbose):
            print("Adding file %s") % filename
        
        if(not os.path.isfile(filename)):
            print("   File %s does not exist!") %filename
            return False

        tdf = TlacDat(filename)

        ## first file added
        if(len(tdl) == 0):
            tdl.append(tdf)
            self.header = tdl[0].header
            return True

        ## check if same run parameters
        if check:
            if(not self.ok_to_add(tdf)):
                if(self.verbose):
                    print "   File *not* added!"
                return False

        tdl.append(tdf)
        
        return True


    def ok_to_add(self, tdf, verbose = False):
        """checks if TlacDat tdf is ok to add to current dat set
        """

        # runtime parameters
        if(not self.header.equal_parameters(tdf.header,
                                            verbose = self.verbose)):
            if(self.verbose):
                print "   Unequal runtime parameters."
            return False

        # random seeds
        rs = []
        for ctd in self.tlac_dat:
            rs.append(ctd.header.get("run", "random_seed"))

        myrs = tdf.header.get("run", "random_seed")
        if(myrs in rs):
            if(self.verbose):
                print "   Random seed %s already there!" %myrs
            return False

        return True

    
    def get_filenames(self):
        """Returns list of tlac dat filenames
        """
        r = []
        for i in self.tlac_dat:
            r.append(i.filename)

        return r
        

    def get_nphot(self, phot_type = "Lya"):
        if(phot_type == "Lya"):
            return np.sum([ i.nphot for i in self.tlac_dat ] )
        elif(phot_type == "UV"):
            return np.sum([ i.nphot_uv for i in self.tlac_dat ] )
        else:
            raise ValueError("Choose \"Lya\" or \"UV\" as `phot_type`!")

    def get_ndat(self, phot_type = "Lya"):
        if(phot_type == "Lya"):
            return np.sum([ i.ndat for i in self.tlac_dat ] )
        elif(phot_type == "UV"):
            return np.sum([ i.ndat_uv for i in self.tlac_dat ] )
        else:
            raise ValueError("Choose \"Lya\" or \"UV\" as `phot_type`!")

        

    def get_EWi(self):
        """Returns intrinsic equivalent width in Angstrom.
        """
        nuv = self.get_nphot('UV')
        if(nuv == 0 ):
            return np.inf

        # Check right frequency modes were used
        assert self.header['emission', 'frequency_mode'] == 2
        assert self.header['emission', 'uv_frequency_mode'] == 3

        deltav =  2. * float(self.header['emission','uv_frequency_param']) * 1e3
        nlya = self.get_nphot('Lya')
        
        return float(nlya) / float(nuv) * deltav / 2.998e8 * 1216.


    def get_v_thermal(self):
        """Returns thermal velocity in m/s.
        """
        return self.header.get_v_thermal()


    def get_frequencies(self, phot_type = 'both', freq_type = 'x'):
        """
        Obtain all frequencies from DatSet as 1D array
        
        Keyword Arguments:
        freq_type -- Defines in what unit the frequencies will be returned
                   * 'x' (default) for x = (nu - nu_0) / Delta nu (code units)
                   * 'v' for velocity shift in km/s
        phot_type -- Can be 'both', 'Lya' or 'UV' depending on what type
        """
        if phot_type == 'both':
            x = np.concatenate((self['Lya','x'], self['UV','x']))
        else:
            x = self[phot_type,'x']
            
        if freq_type == 'x':
            return x
        elif freq_type == 'v':
            v_th = self.get_v_thermal()
            c = 2.998e8
            return -c / (1. + c / (v_th * x)) / 1e3
        else:
            raise ValueError("freq_type can be 'x' or 'v'.")
            
            
               
    def load(self, phot_type = 'both', *args):
        for cf in self.tlac_dat:
            cf.load(phot_type, *args)
    
    
    def load_all(self, phot_type = 'both', **kwargs):
        """Loads all photons of `phot_type`.
        phot_type can be 'both', 'Lya' or 'UV'.
        """
        for cf in self.tlac_dat:
            cf.load_all(phot_type, **kwargs)


    def unload_all(self, **kwargs):
        for cf in self.tlac_dat:
            cf.unload_all(**kwargs)
        


    
