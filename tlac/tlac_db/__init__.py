
import os, sys, inspect
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(
    inspect.getfile( inspect.currentframe() ))[0],"lib")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
del cmd_subfolder

from TlacDB import TlacDB
from TlacDatSet import TlacDatSet
from TlacHeader import TlacHeader
from TlacDatFile import TlacDatFile
