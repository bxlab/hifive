import sys

import _bi
import _fivec_binning
import _fivec_optimize
import _hic_binning
import _hic_distance
import _hic_interactions
import _hic_optimize


__version__ = '2.0.0'
__author__ = 'Michael Sauria'

#------------------------------------------------

class Std( object ):
    def __init__( self ):
        self.__dict__ = sys.modules
    def __getattr__( self, name ):
        m = __import__(name)
        return m

std = Std()
