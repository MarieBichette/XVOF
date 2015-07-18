#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
'''
Test de performance sur le calcul en formulation VE de l'équation d'état de type Mie-Gruneisen
'''


if __name__ == '__main__':
    #
    # Lancement du profiling
    #
    from xvof.utilities import timeit_file
    from os import system
    from random import random
    from xvof.equationsofstate.miegruneisen import MieGruneisen
    #

    @timeit_file('/tmp/solve_ve.log')
    def profil_solve_ve(equation_of_state, spec_vol, e_int, it_nb=100):
        """
        Fait it_nb appel(s) à  solve_ve à des fins de profiling
        """
        i = 0
        while i < it_nb:
            equation_of_state.solveVolumeEnergy(spec_vol, e_int)
            i += 1
    #
    MY_EE = MieGruneisen()
    RHO = 8000.0 + random() * 1000.0
    E_INT = random() * 2.0e+03
    NBR_ITER = 1110000
    print "Lancement profiling sur {:6d} itÃ©rations".format(NBR_ITER)
    profil_solve_ve(MY_EE, 1. / RHO, E_INT, NBR_ITER)
    system('cat /tmp/solve_ve.log')
