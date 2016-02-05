import numpy as np

from TlacFile import TlacFile

class TlacASCII(TlacFile):
    """For reading tlac *dat and *uvdat files."""

    def __init__(self, filename):
        TlacFile.__init__(self, filename)
        self.filename = filename
        
        self.rawdat = []
        self.dat = {}
        
        self.icollst = []
        self.set_icollst()
        
    def __getitem__(self, key):
        if(not key in self.dat):
            self.load(key)
        return self.dat[key]

    def keys(self):
        return self.icollst
        
    def load(self, keys):
        if(not isinstance(keys, list)):
            keys = [keys ]

        icols = self.keys_to_icols(keys)
        if(len(icols) == len(self.icollst)):
            icols = None # load all
            
        self.rawdat.append(np.loadtxt(self.filename, usecols = icols))

        for i,k in enumerate(keys):
            if(k in self.dat):
                raise Exception("Key " + k + " already loaded")
            if(len(keys) == 1):
                self.dat[k] = self.rawdat[-1]
            else:
                self.dat[k] = self.rawdat[-1][:,i]

    
    def load_all(self):
        self.load(self.icollst)
        
    def keys_to_icols(self, keys):        
        """Given a list of keys, returns in which columns they are in the file.
        """
        return [ self.icollst.index(k) for k in keys ]
        
    def set_icollst(self):
        """Sets column headers according to file's first line.
        """
        f = open(self.filename, "r")
        l = f.readline()
        f.close()

        if(l[:1] != "#"):
            raise Exception("Expected first line commented in " + self.filename)
            
        l = l[1:].strip()
        cols = l.split(" ")

        r = []
        for i, ccol in enumerate(cols):
            name, rest = ccol.split("(")
            name = name.strip()
            rest = rest[:-1].split("-")
            if(len(rest) > 1):
                a = int(rest[0])
                b = int(rest[1])
                for j in range(0, b - a + 1):
                    r.append(name + str(j))
            else:
                r.append(name)


        self.icollst = r
