import logging
from fastlmm.util.runner import *
from contextlib import contextmanager
import threading

dyn = threading.local()

# from short example in http://stackoverflow.com/questions/2001138/how-to-create-dynamical-scoped-variables-in-python999
@contextmanager
def _dyn_vars(**new):
    old = {}
    for name, value in new.items():
        old[name] = getattr(dyn, name, None)
        setattr(dyn, name, value)
    yield
    for name, value in old.items():
        setattr(dyn, name, value)

class _MapReduce(object): #implements IDistributable
    """
    class to run distributed map using the idistributable back-end
    """


    def __init__(self, input_seq, mapper, nested, reducer, input_files=None, output_files=None, name=None):

        self.input_seq = input_seq
        self.mapper = mapper
        self.nested = nested
        if (self.mapper is not _identity) and (self.nested is not None):
            raise Exception("'mapper' and 'nested' should not both be set")
        self.reducer = reducer
        self.name = name

        if input_files is None:
            self.input_files = []
        else:
            self.input_files = input_files

        if output_files is None:
            self.output_files = []
        else:
            self.output_files = output_files


#start of IDistributable interface--------------------------------------
    @property
    def work_count(self):
        return len(self.input_seq)

    def work_sequence_range(self, start, stop):
        for i in xrange(start,stop):
            input_arg = self.input_seq[i]
            if self.nested is None:
                #logging.debug("\nrandom access executing %i" % i)
                with _dyn_vars(is_in_nested=False):
                    yield lambda i=i, input_arg=input_arg: self.dowork(i, input_arg)   # the 'i=i',etc is need to get around a strangeness in Python
            else:
                assert self.nested is not None, "real assert"
                with _dyn_vars(is_in_nested=True):
                    dist = self.nested(input_arg)
                    yield dist

    def work_sequence(self):
        for i, input_arg in enumerate(self.input_seq):
            if self.nested is None:
                #logging.debug("\nexecuting %i" % i)
                with _dyn_vars(is_in_nested=False):
                    yield lambda i=i, input_arg=input_arg: self.dowork(i, input_arg)  # the 'i=i',etc is need to get around a strangeness in Python
            else:
                assert self.nested is not None, "real assert"
                with _dyn_vars(is_in_nested=True):
                    dist = self.nested(input_arg)
                    yield dist


    def reduce(self, output_seq):
        '''
        '''

        return self.reducer(output_seq)


    #optional override -- the str name of the instance is used by the cluster as the job name
    def __str__(self):
        if self.name is None:
            return "map_reduce()"
        else:
            return self.name
 #end of IDistributable interface---------------------------------------

    def dowork(self, i, input_arg):
        #logging.info("{0}, {1}".format(len(train_snp_idx), len(test_snp_idx)))
        #logging.debug("\nexecuting {0}".format(input_arg))
        work = lambda : self.mapper(input_arg)
        result = run_all_in_memory(work)
        return result

   
    # required by IDistributable
    @property
    def tempdirectory(self):
        return ".work_directory.{0}".format(self.name)
        

    def copyinputs(self, copier):
        for fn in self.input_files:
            copier.input(fn)

    def copyoutputs(self,copier):
        for fn in self.output_files:
            copier.output(fn)

def _identity(x):
    return x

def _is_in_nested():
    return hasattr(dyn,"is_in_nested") and dyn.is_in_nested

def map_reduce(input_seq, mapper=_identity, reducer=list, input_files=None, output_files=None, name=None, runner=None, nested=None):
    """
    Function for running a function on sequence of inputs and running a second function on the results. Can be nested and clusterized.

    :param input_seq: a sequence of inputs. The sequence must support the len function and be indexable. e.g. a list, xrange(100)
    :type input_seq: a sequence

    :param mapper: A function to apply to each set of inputs (optional). Defaults to the identity function. (Also see 'mapper')
    :type mapper: a function

    :param reducer: A function to turn the results from the mapper to a single value (optional). Defaults to creating a list of the results.
    :type reducer: a function that takes a sequence

    :param input_files: A list that tells what input files are needed. The list can contain the names of files (strings), None (ignored), or
        objects such as :class:`.SnpReader`'s that can self-report their input files.
    :type input_files: a list

    :param output_files: A list that tells what output files will be produced. The list can contain the names of files (strings), None (ignored), or
        objects such as :meth:`.map_reduce`'s that can self-report their output files.
    :type output_files: a list

    :param name: A name to be displayed if this work is done on a cluster.
    :type name: a string

    :param runner: a runner, optional: Tells how to run locally, multi-processor, or on a cluster.
        If not given, the function is run locally.
    :type runner: a runner.

    :param nested: a mapper function that is itself a map_reduce (or other IDistributable). Some runners can efficiently clusterize such nested mappers. 
    :type nested: a function

    :rtype: The results from the reducer.

    :Example:

    Square the numbers 1 to n (inclusive) and report their sum.

        >>> def sos(n):
        ...   return map_reduce(xrange(1,n+1), 
        ...                     mapper=lambda x: x*x,
        ...                     reducer=sum)
        >>> sos(99)
        328350

    Make a list of sos(0) to sos(9). By using nested instead of mapper, we telling runners they can try to distribute both levels of map_reduce.

        >>> map_reduce(xrange(10),
        ...            nested=lambda s:sos(s),
        ...            reducer=list)
        [0, 1, 5, 14, 30, 55, 91, 140, 204, 285]

    """

    dist = _MapReduce(input_seq, mapper=mapper, nested=nested, reducer=reducer, input_files=input_files, output_files=output_files,name=name)
    if runner is None and _is_in_nested():
        return dist

    if runner is None:
        runner = Local()

    result = runner.run(dist)
    return result
    
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
