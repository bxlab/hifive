import sys

import binning
import data
import _binning
import _distance


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
