#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-

from os import system

from xvof.element.element1d import Element1d
from xvof.equationsofstate.miegruneisen import MieGruneisen
from xvof.utilities import timeit_file


if __name__ == '__main__':
    #
    # Lancement du profiling
    #
    @timeit_file("/tmp/newton_raphson.profil")
    def profil_newton_raphson(eos, rho_old, rho_new, p_old, pseudo_old, nrj_old, nit=100):
        it = 0
        while(it < nit):
            Element1d.newton_raphson_for_ve(eos, rho_old, rho_new, p_old, pseudo_old, nrj_old)
            it += 1

    EE = MieGruneisen()
    NBR_ITER = 375000
    print "Lancement profiling sur {:6d} itérations".format(NBR_ITER)
    profil_newton_raphson(EE, 8500., 8600., 1.0e+09, 0., 1.e+04, NBR_ITER)
    system('cat /tmp/newton_raphson.profil')