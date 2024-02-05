import sys

from . import _fivec_binning
from . import _fivec_optimize
from . import _hic_binning
from . import _hic_distance
from . import _hic_interactions
from . import _hic_optimize
from . import _hic_domains
from . import _hmm
from . import _quasar
from . import hmm

#------------------------------------------------

class Std( object ):
    def __init__( self ):
        self.__dict__ = sys.modules
    def __getattr__( self, name ):
        m = __import__(name)
        return m

std = Std()
