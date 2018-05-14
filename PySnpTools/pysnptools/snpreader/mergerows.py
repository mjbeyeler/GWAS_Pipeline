#import numpy as np
#import subprocess, sys, os.path
#from itertools import *
#import pandas as pd
#import logging
#from snpreader import SnpReader
#from pysnptools.standardizer import Unit
#from pysnptools.standardizer import Identity
#from pysnptools.pstreader import PstData
#import warnings
#import time

#def _iidset(reader):
#    if reader is None:
#        return set()
#    return {tuple(item) for item in reader.iid}

#def _in_other(i_little,iidset_list):
#    little = iidset_list[i_little]
#    for i_big,big in enumerate(iidset_list):
#        #Set 'a' is in set 'b' if
#        #         'a' is a proper subset of 'b'
#        #         'a' equals 'b' but 'b' is listed first (this also stops a set from being a subset of itself)
#        if little < big or (i_big < i_little and little == big):
#            return True
#    return False

##!!!this is unused and untested. Also, since it only works with SnpReader a better name would be MergeByIid
## should we first confirm that all col_poperty values match across the items? (if so, do NaN, right too)
#class _MergeRows(SnpReader):
#    @staticmethod
#    def factory(*readerlist):
#        #Remove any readers for which another reader has all the same row ids
#        iidset_list = [_iidset(reader) for reader in readerlist]
#        readerlist = [reader for index, reader in enumerate(readerlist) if not _in_other(index,iidset_list)]

#        if len(readerlist) == 0:
#            return None
#        if len(readerlist) == 1:
#            return readerlist[0]
        
#        return _MergeRows(*readerlist)

#    def __init__(self, *readerlist):
#        self.readerlist = readerlist

#    def __repr__(self): 
#        return "{0}({1})".format(self.__class__.__name__,",".join([repr(reader) for reader in self.readerlist]))

#    @property
#    def row(self):
#        if not hasattr(self,"_row"):
#            self._row = np.concatenate([reader.row for reader in self.readerlist])
#        return self._row

#    @property
#    def col(self):
#        if not hasattr(self,"_col"):
#            self._col = self.readerlist[0].col
#            for i in xrange(1,len(self.readerlist)):
#                assert np.array_equal(self._col,self.readerlist[i].col), "all col's must be the same"
#        return self._col

#    @property
#    def col_property(self):
#        return self.readerlist[0].col_property

#    def _find_one(self,iid_index_or_none,sid_index_or_none):
#        assert sid_index_or_none is None, "Expect sid_index_or_none to be None"
#        assert iid_index_or_none is not None, "Expect iid_index_or_none to be not None"

#        result = None
#        iid_goal = self.iid[iid_index_or_none]
#        for i, reader in enumerate(self.readerlist):
#            try:
#                iididx = reader.iid_to_index(iid_goal)
#            except:
#                continue # leave the loop
#            assert result is None or len(result[1])==0, "for now code assumes all values will be read from one part of merged SnpReader"
#            result = i, iididx
#        assert result is not None and len(result[1]) == len(iid_goal), "Could not find all indexes."
#        return result


#    def _read(self, iid_index_or_none, sid_index_or_none, order, dtype, force_python_only, view_ok):
#        i, iididx = self._find_one(iid_index_or_none, sid_index_or_none)
#        result = readerlist[i]._read(iididx, sid_index_or_none, order, dtype, force_python_only, view_ok)
#        return result

#    def __getitem__(self, iid_indexer_and_snp_indexer):
#        if isinstance(iid_indexer_and_snp_indexer,tuple): # similar code elsewhere
#            iid0_indexer, iid1_indexer = iid_indexer_and_snp_indexer
#        else:
#            iid0_indexer = iid_indexer_and_snp_indexer
#            iid1_indexer = iid0_indexer

#        i, iididx = self._find_one(iid0_indexer, None)
#        result = self.readerlist[i][iididx,iid1_indexer]
#        return result
