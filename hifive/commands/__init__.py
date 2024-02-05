import sys

from . import connect_files
from . import create_fragmentset
from . import create_fivec_dataset
from . import create_fivec_project
from . import normalize_fivec_project
from . import complete_fivec_project
from . import create_fivec_heatmap
from . import combine_fivec_replicates
from . import get_fivec_interval
from . import create_fendset
from . import create_hic_dataset
from . import create_hic_project
from . import normalize_hic_project
from . import complete_hic_project
from . import create_hic_heatmap
from . import create_hic_mrheatmap
from . import combine_hic_replicates
from . import get_hic_interval
from . import find_quasar_scores

#------------------------------------------------

class Std( object ):
    def __init__( self ):
        self.__dict__ = sys.modules
    def __getattr__( self, name ):
        m = __import__(name)
        return m

std = Std()
