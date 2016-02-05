from TlacFile import TlacFile
import re

class TlacDatFile(TlacFile):
    def __init__(self, filename):
        """Loads *.dat files produced by tlac. Loads header of file."""
        TlacFile.__init__(self, filename)

        

    def add_data(self, tlac_file):
        TlacFile.add_data(self, tlac_file)

        self.Nphot += tlac_file.Nphot
        self.fesc = self.ndat / self.Nphot
        
        
