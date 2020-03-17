#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Module simpliste pour du profiling simpliste
"""

import os
import sys
from collections import OrderedDict
from time import clock, time, sleep

cumul_times = OrderedDict()


def timeit_file(filename=None):
    """
    Decorateur permettant de sortir dans le fichier de nom "filename"
    le temps d'éxécution réel et cpu de la fonction décorée

    @param
    filename : nom du fichier de sortie
    """
    def timeit(func):
        """
        Renvoi la fonction wrapper
        """
        def wrapper(*args, **kwargs):
            """
            wrapper autour de la fonction func
            """
            if filename is not None:
                file_desc = open(filename, 'a+')
            else:
                file_desc = sys.stdout
            begin_cpu_time = clock()
            begin_real_time = time()
            print >> file_desc, "=" * 80
            print >> file_desc, "Appel de {:s}".format(func.__name__)
            fonction = func(*args, **kwargs)
            cpu_duration = clock() - begin_cpu_time
            real_duration = time() - begin_real_time
            record = cumul_times.setdefault(func.__name__, [0., 0.])
            record[0] += cpu_duration
            record[1] += real_duration
            print >> file_desc, "\t Durée CPU instantanée/cumulée : {:4f}/{:4f}s".format(cpu_duration, record[0])
            print >> file_desc, "\t Durée réelle instantanée/cumulée : {:4f}/{:4f}s".format(real_duration, record[1])
            print >> file_desc, "\t Durée CPU/réelle totale : {:4f}/{:4f}" \
                .format(sum([rec[0] for rec in cumul_times.values()]), sum([rec[1] for rec in cumul_times.values()]))
            print >> file_desc, "=" * 80 + os.linesep
            if filename is not None:
                file_desc.close()
            return fonction
        return wrapper
    return timeit


def logit(func):
    """
    Décorateur permettant de tracer l'appel des fonctions décorées en affichant
    leurs argurments

    @param func : fonction dont on veut tracer l'appel
    """
    def wrapper(*args, **kwargs):
        """
        Wrapper autour de la fonction func
        """
        print "Appel de {:s} ".format(func.__name__,)
        print "Arguments : {:s}".format(args)
        print "Arguments par mots clefs : {:s}".format(kwargs)
        fonction = func(*args, **kwargs)
        return fonction
    return wrapper


if __name__ == "__main__":
    @timeit_file('toto.log')
    @logit
    def countdown(nbr_iter):
        """
        Compte à rebours sur nbr_iter itérations
        """
        for i in range(nbr_iter, 0, -1):
            if i % 1000 == 0:
                print i
            sleep(0.0001)

    countdown(10000)
    countdown(10000)
    print "name : ", countdown.__name__
