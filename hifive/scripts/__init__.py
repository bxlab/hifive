import sys

from . import fetch_mrh_data
from . import hifive
from . import hifive2butlr
from . import hifive2cooler
from . import hifive2mcool

#------------------------------------------------

class Std( object ):
    def __init__( self ):
        self.__dict__ = sys.modules
    def __getattr__( self, name ):
        m = __import__(name)
        return m

std = Std()
