#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Module simpliste pour du profiling simpliste
"""

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############ IMPORTATIONS DIVERSES  ####################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
from time import clock, time, sleep
import sys

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
####### DEFINITION DES CLASSES & FONCTIONS  ###############
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


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
            file_desc = None
            if filename is not None:
                file_desc = open(filename, 'a+')
            else:
                file_desc = sys.stdout
            begin_cpu_time = clock()
            begin_real_time = time()
            print >> file_desc, "=" * 80
            print >> file_desc, "Appel de {:s}".format(func.__name__)
            fonction = func(*args, **kwargs)
            end_cpu_time = clock()
            end_real_time = time()
            print >> file_desc, "\t Durée CPU : {:4f}s".\
            format(end_cpu_time - begin_cpu_time)
            print >> file_desc, "\t Durée réelle : {:4f}s".\
            format(end_real_time - begin_real_time)
            print >> file_desc, "=" * 80 + "\n"
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############ PROGRAMME PRINCIPAL ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if __name__ == "__main__":
    @timeit_file('toto.log')
    @logit
    def countdown(nbr_iter):
        """
        Compte à rebours sur nbr_iter itérations
        """
        for i in range(nbr_iter, 0, -1):
            if (i % 1000 == 0):
                print i
            sleep(0.0001)

    countdown(10000)
    countdown(10000)
    print "name : ", countdown.__name__