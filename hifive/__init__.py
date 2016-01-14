import sys

import fend
import fragment
import fivec_data
import hic_data
import fivec
import hic
import fivec_binning
import hic_binning
import plotting
import hic_domains

from hic import HiC
from fivec import FiveC
from fend import Fend
from fragment import Fragment
from hic_data import HiCData
from fivec_data import FiveCData
from hic_domains import TAD
from hic_domains import Compartment

from .version import version as __version__

#------------------------------------------------

class Std( object ):
    def __init__( self ):
        self.__dict__ = sys.modules
    def __getattr__( self, name ):
        m = __import__(name)
        return m

std = Std()
