import os.path
import numpy as np

class TlacFile:
    def __init__(self, filename):
        """The class is not to be importded directly.
        It's a general class where specific TlacTYPE can inherit from.
        """
        self.filename = os.path.expanduser(filename)

        if(not os.path.isfile(self.filename)):
            raise Exception("File " + self.filename + " does not exist!")


        
