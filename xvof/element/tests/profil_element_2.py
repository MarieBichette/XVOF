#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-

import numpy as np
from os import system

from xvof.element.one_dimension_element import OneDimensionElement
from xvof.equationsofstate.miegruneisen import MieGruneisen
from xvof.miscellaneous import *
from xvof.node.node1d import Node1d
from xvof.utilities import timeit_file

if __name__ == '__main__':
    #
    # Lancement du profiling
    #
    @timeit_file("/tmp/nouvo_pression.profil")
    def profil_calculer_nouvo_pression(nit=100):
        it = 0
        equation_detat = MieGruneisen()
        num_props = numerical_props(0.2, 1.0, 0.35)
        mat_props = material_props(1.0e+05, 0.0, 8129., equation_detat)
        geom_props = geometrical_props(1.0e-06)
        props = properties(num_props, mat_props, geom_props)
        noda = Node1d(1, np.array([4.5e-03]))
        nodb = Node1d(2, np.array([7.0e-03]))
        noda._xtpdt = np.array([5.0e-03])
        nodb._xtpdt = np.array([6.25e-03])
        my_elem = OneDimensionElement(props, 1, [noda, nodb])
        while(it < nit):
            my_elem.calculer_nouvo_pression()
            it += 1

    EE = MieGruneisen()
    NBR_ITER = 375000
    print "Lancement profiling sur {:6d} itérations".format(NBR_ITER)
    profil_calculer_nouvo_pression(NBR_ITER)
    system('cat /tmp/nouvo_pression.profil')