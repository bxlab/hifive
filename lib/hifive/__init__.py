import sys

import hic
import fivec
import fend
import fragment
import analysis
import plotting
import bi
import _bi


"""
Shamelessly ripped off of py.std
"""
__version__ = '0.0.0'
__author__ = 'Michael Sauria'

#------------------------------------------------

class Std( object ):
    def __init__( self ):
        self.__dict__ = sys.modules
    def __getattr__( self, name ):
        m = __import__(name)
        return m

std = Std()
