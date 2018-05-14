import numpy as np
import scipy as sp
import logging
from pysnptools.standardizer import Standardizer
import warnings

class BySidCount(Standardizer):
    '''
    '''

    def __init__(self):
        super(BySidCount, self).__init__()
        warnings.warn("BySidCount no longer supported", DeprecationWarning)

