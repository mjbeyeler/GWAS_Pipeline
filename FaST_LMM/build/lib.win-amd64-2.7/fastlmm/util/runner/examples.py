import math
from fastlmm.util.mapreduce import map_reduce
from fastlmm.util.runner import Local, LocalMultiProc, HPC, LocalMultiThread
import os

def is_prime(n):
    assert n == int(n) and n>1, "Expect integers greater than 1"
    for j in xrange(2,int(math.sqrt(n))+1):
        if n % j == 0:
            return False
    return True

# Iterative algorithm for finding prime numbers in a range
def prime_search0(start,stop):
    assert start < stop, "start must be less than stop"
    prime_list = []
    for i in xrange(start,stop):
        if is_prime(i):
            prime_list.append(i)
    return prime_list

# The similar map_reduce algorithm for finding prime numbers in a range
def prime_search1(start,stop,runner):

    def mapper(i):
        if is_prime(i):
            #assert i != 5, "I just don't like fives"
            return i
        else:
            return None

    def reducer(sequence):
        result = []
        for i in sequence:
            if i is not None:
                result.append(i)
        return result

    return map_reduce(xrange(start,stop),
                        mapper=mapper,
                        reducer=reducer, #lambda sequence: [i for i in sequence if i is not None], #Filter out the None's
                        runner=runner)

if __name__ == '__main__':
    #Run the iterative algorithm
    #print prime_search0(2,10) #=> [2, 3, 5, 7]

    #Run the map_reduce algorithm locally.
    #print prime_search1(2,10,runner=Local()) #=> [2, 3, 5, 7]

    #Now we run map_reduce on 20 processors.
    #from PrimeNumbers.examples import prime_search1 #If not running local, must import your function. (Recall that for import to work, you also need an empty __init__.py).
    #print prime_search1(2,10,runner=LocalMultiProc(20)) #=> [2, 3, 5, 7]

    #Finally we run on HPC
    #------- To run on HPC must create an hpc cluster object
    #remote_python_parent=r"\\GCR\Scratch\B99\escience\{0}\ppv0".format(os.environ['USERNAME']) #where to copy your "python_path" code to.
    #hpc_runner= HPC(10, 'GCR',r"\\GCR\Scratch\B99\escience",
    #                                        remote_python_parent=remote_python_parent,
    #                                        unit='node', #core, socket, node
    #                                        update_remote_python_parent=True,
    #                                        template="Preemptable",
    #                                        priority="Lowest",
    #                                        nodegroups="Preemptable",
    #                                        runtime="0:11:0", # day:hour:min
    #                                        )
    #runner=LocalMultiProc(2,just_one_process=False)
    #runner = Local()
    runner = LocalMultiThread(2,just_one_process=False)
    print prime_search1(2,10,runner=runner) #=> [2, 3, 5, 7]
    print "done"