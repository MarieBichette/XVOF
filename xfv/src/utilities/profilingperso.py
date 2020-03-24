#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A simple module for simple profiling
"""

import os
import sys
from collections import OrderedDict
from time import perf_counter, time, sleep


CUMUL_TIMES = OrderedDict()


def timeit_file(filename=None):
    """
    A decorator allowing to write in the file, the real and cpu execution
    times of the decorated function

    @param
    filename : name of the output file
    """
    def timeit(func):
        """
        Returns the wrapper function
        """
        def wrapper(*args, **kwargs):
            """
            Wrapper around the function func
            """
            if filename is not None:
                file_desc = open(filename, 'a+')
            else:
                file_desc = sys.stdout
            begin_cpu_time = perf_counter()
            begin_real_time = time()
            print("=" * 80, file=file_desc)
            print("Call of {:s}".format(func.__name__), file=file_desc)
            fonction = func(*args, **kwargs)
            cpu_duration = perf_counter() - begin_cpu_time
            real_duration = time() - begin_real_time
            record = CUMUL_TIMES.setdefault(func.__name__, [0., 0.])
            record[0] += cpu_duration
            record[1] += real_duration
            print("\t Instantaneous/cumulative CPU duration: {:4f}/{:4f}s"
                  .format(cpu_duration, record[0]), file=file_desc)
            print("\t Instantaneous/cumulative real duration: {:4f}/{:4f}s"
                  .format(real_duration, record[1]), file=file_desc)
            print("\t Total real CPU duration: {:4f}/{:4f}"
                  .format(sum([rec[0] for rec in list(CUMUL_TIMES.values())]),
                          sum([rec[1] for rec in list(CUMUL_TIMES.values())])),
                  file=file_desc)
            print("=" * 80 + os.linesep, file=file_desc)
            if filename is not None:
                file_desc.close()
            return fonction
        return wrapper
    return timeit


def logit(func):
    """
    A decorator allowing to trace the call of decorated function by printing
    their arguments

    @param func : function which call should be traced
    """
    def wrapper(*args, **kwargs):
        """
        Wrapper around the function func
        """
        print("Call of {:s} ".format(func.__name__,))
        print("Positional arguments: {}".format(args))
        print("Keyword arguments: {}".format(kwargs))
        fonction = func(*args, **kwargs)
        return fonction
    return wrapper


if __name__ == "__main__":
    @timeit_file('toto.log')
    @logit
    def countdown(nbr_iter):
        """
        Countdown on nbr_iter iterations
        """
        for i in range(nbr_iter, 0, -1):
            if i % 1000 == 0:
                print(i)
            sleep(0.0001)

    countdown(10000)
    countdown(10000)
    print("name : ", countdown.__name__)
