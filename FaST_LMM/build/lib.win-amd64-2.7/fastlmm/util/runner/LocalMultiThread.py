'''
Runs a distributable job on multiple processors.  Returns the value of the job.

See SamplePi.py for examples.
'''

from fastlmm.util.runner import *
import os
import logging
try:
    import dill as pickle
except:
    logging.warning("Can't import dill, so won't be able to clusterize lambda expressions. If you try, you'll get this error 'Can't pickle <type 'function'>: attribute lookup __builtin__.function failed'")
    import cPickle as pickle
import subprocess, sys, os.path
import threading
import fastlmm.util.util as util
from Queue import PriorityQueue

class LocalMultiThread: # implements IRunner
    '''Designed so that reduce will start running as soon as the 1st task as finished
    '''

    def __init__(self, taskcount, mkl_num_threads = None, just_one_process = False,):
        if not 0 < taskcount: raise Exception("Expect taskcount to be positive")

        self.taskcount = taskcount
        self.just_one_process = just_one_process
        if mkl_num_threads != None:
            os.environ['MKL_NUM_THREADS'] = str(mkl_num_threads)

    def _result_sequence(self,thread_list,priority_queue,shaped_distributable):
        for thread in thread_list:
            if not self.just_one_process:
                thread.join()
            result_sequence = priority_queue.get()[1]
            for result in result_sequence:
                yield result

    def run(self, distributable):
        JustCheckExists().input(distributable)

        priority_queue = PriorityQueue()
        thread_list = []
        shaped_distributable = shape_to_desired_workcount(distributable, self.taskcount)
        for taskindex in xrange(self.taskcount):
            def _target(taskindex=taskindex):
                result_list = []
                for work in work_sequence_for_one_index(shaped_distributable, self.taskcount, taskindex):
                    result_list.append(run_all_in_memory(work))
                priority_queue.put((taskindex,result_list))
            if not self.just_one_process:
                thread = threading.Thread(target=_target,name=str(taskindex))
                thread_list.append(thread)
                thread.start()
            else:
                thread_list.append(None)
                _target()
        
        result_sequence = self._result_sequence(thread_list, priority_queue,shaped_distributable)
        result = shaped_distributable.reduce(result_sequence)

        JustCheckExists().output(distributable)
        return result
