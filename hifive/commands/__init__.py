import sys

import connect_files
import create_fragmentset
import create_fivec_dataset
import create_fivec_project
import normalize_fivec_project
import complete_fivec_project
import create_fivec_heatmap
import combine_fivec_replicates
import get_fivec_interval
import create_fendset
import create_hic_dataset
import create_hic_project
import normalize_hic_project
import complete_hic_project
import create_hic_heatmap
import create_hic_mrheatmap
import combine_hic_replicates
import get_hic_interval

#------------------------------------------------

class Std( object ):
    def __init__( self ):
        self.__dict__ = sys.modules
    def __getattr__( self, name ):
        m = __import__(name)
        return m

std = Std()
