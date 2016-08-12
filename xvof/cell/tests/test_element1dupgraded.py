#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module element1dupgraded
"""
import numpy as np
import unittest

import xvof.cell.element1dupgraded as el1dup
from xvof.miscellaneous import geometrical_props, material_props
from xvof.miscellaneous import numerical_props, properties

import xvof.cell.one_dimension_cell as el1d
from xvof.equationsofstate import MieGruneisen
from xvof.node import Node1d


class Element1dUpgradedTest(unittest.TestCase):

    def setUp(self):
        """ Préparation des tests """
        ee = MieGruneisen()
        num_props = numerical_props(0.2, 1.0, 0.35)
        mat_props = material_props(1.0e+05, 0.0, 8129., ee)
        geom_props = geometrical_props(1.0e-06)
        props = properties(num_props, mat_props, geom_props)
        noda = Node1d(1, poz_init=np.array([0.6]), section=1.0e-06)
        nodb = Node1d(1, poz_init=np.array([-0.2]), section=1.0e-06)
        my_elem = el1d.OneDimensionCell(props, 123, [noda, nodb])
        self.my_elem_up = el1dup.Element1dUpgraded(my_elem, 0.5)

    def tearDown(self):
        pass

    def test_calculer_nouvo_taille(self):
        """ Test de la méthode Element1dUpgraded.calculer_nouvo_taille """
        self.my_elem_up.calculer_nouvo_taille(1.0e-06)
        self.assertEqual(self.my_elem_up.taille_t_plus_dt_gauche, 0.4)
        self.assertEqual(self.my_elem_up.taille_t_plus_dt_droite, 0.4)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
