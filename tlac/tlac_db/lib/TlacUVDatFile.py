from TlacFile import TlacFile
import re

class TlacUVDatFile(TlacFile):
    def __init__(self, filename):
        """Loads *.uvdat files produced by tlac.
        """
        TlacFile.__init__(self, filename)

        self.nphot = float(self.header.get("run", "nphot_uv"))
        self.ndat = float(self.header.get("result", "ndat_uv"))

        
