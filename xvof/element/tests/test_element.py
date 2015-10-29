#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module element
"""
import unittest
import numpy as np

import xvof.element.element as el
from xvof.equationsofstate import MieGruneisen
from xvof.miscellaneous import geometrical_props, material_props
from xvof.miscellaneous import numerical_props, properties
from xvof.node import Node1d


class ElementTest(unittest.TestCase):

    def setUp(self):
        """ Préparation des tests """
        equation_detat = MieGruneisen()
        num_props = numerical_props(0.2, 1.0, 0.35)
        mat_props = material_props(1.0e+05, 0.0, 8129., equation_detat)
        geom_props = geometrical_props(1.0e-06)
        props = properties(num_props, mat_props, geom_props)
        noda = Node1d(1, np.array([-0.5]))
        nodb = Node1d(2, np.array([0.1]))
        self.my_elem = el.Element(props, 1, [noda, nodb])

    def tearDown(self):
        pass

    def test_coord(self):
        """ Test de la méthode Element.coord() """
        np.testing.assert_array_equal(self.my_elem.coord, np.array([-0.2]))

    def test_incrementer(self):
        """ Test de la méthode Element.incrementer() """
        self.my_elem._pression_t_plus_dt = 2.0e+09
        self.my_elem._rho_t_plus_dt = 8500.
        self.my_elem._cson_t_plus_dt = 4360.
        self.my_elem._nrj_t_plus_dt = 1.0e+06
        self.my_elem._size_t_plus_dt = 0.015
        self.my_elem.incrementer()
        self.assertEqual(self.my_elem.pression_t, 2.0e+09)
        self.assertEqual(self.my_elem.rho_t, 8500.)
        self.assertEqual(self.my_elem.cson_t, 4360.)
        self.assertEqual(self.my_elem.nrj_t, 1.0e+06)
        self.assertEqual(self.my_elem.taille_t, 0.015)

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
