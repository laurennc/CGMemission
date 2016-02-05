"""
TlacAnalysis

Library to analyze and plot tlac output.
Relies on TlacDB for sorting of data.
"""

# Import TlacDB 
try:
    from TlacDB import TlacDB
except:
    import os
    import sys
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                 '../tlac_db/lib/')))
    from TlacDB import TlacDB


import weighted_kde

from .tlac_plot import *
from .tlac_weights import *
from .spectra import *

from .physics import *
from .analytic_solutions import *

from .grid import *



