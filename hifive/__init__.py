import sys

from . import fend
from . import fragment
from . import fivec_data
from . import hic_data
from . import fivec
from . import hic
from . import fivec_binning
from . import hic_binning
from . import plotting
from . import hic_domains
from . import quasar
from . import schic
from . import scripts
from . import commands
from . import version

from .hic import HiC
from .fivec import FiveC
from .fend import Fend
from .fragment import Fragment
from .hic_data import HiCData
from .fivec_data import FiveCData
from .hic_domains import TAD
from .hic_domains import Compartment
from .quasar import Quasar
from .schic import scHiC

from .version import __version__

#------------------------------------------------

class Std( object ):
    def __init__( self ):
        self.__dict__ = sys.modules
    def __getattr__( self, name ):
        m = __import__(name)
        return m

std = Std()
