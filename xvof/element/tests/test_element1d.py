#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module element1d
"""
import numpy as np
import unittest

import xvof.element.one_dimension_element as el1d
from xvof.equationsofstate import MieGruneisen
from xvof.miscellaneous import geometrical_props, material_props
from xvof.miscellaneous import numerical_props, properties
from xvof.node import Node1d


class Element1dTest(unittest.TestCase):

    def setUp(self):
        """ Préparation des tests """
        equation_detat = MieGruneisen()
        num_props = numerical_props(0.2, 1.0, 0.35)
        mat_props = material_props(1.0e+05, 0.0, 8129., equation_detat)
        geom_props = geometrical_props(1.0e-06)
        props = properties(num_props, mat_props, geom_props)
        noda = Node1d(1, np.array([4.5e-03]))
        nodb = Node1d(2, np.array([7.0e-03]))
        noda._xtpdt = np.array([5.0e-03])
        nodb._xtpdt = np.array([6.25e-03])
        self.my_elem = el1d.OneDimensionElement(props, 1, [noda, nodb])

    def tearDown(self):
        pass

    def test_calculer_nouvo_pression(self):
        """ Test de la méthode Element1d.calculer_nouvo_pression() """
        self.my_elem._rho_t_plus_dt = 9000.0
        self.my_elem.calculer_nouvo_pression()
        self.assertEqual(self.my_elem.nrj_t_plus_dt / 1.e+05,
                         1.0337842440399707)
        self.assertEqual(self.my_elem.pression_t_plus_dt / 1.e+09,
                         17.366763163767697)
        self.assertEqual(self.my_elem.cson_t_plus_dt / 1.0e+03,
                         4.89803134404931)

    def test_calculer_nouvo_taille(self):
        """ Test de la méthode Element1d.calculer_nouvo_taille() """
        self.my_elem.calculer_nouvo_taille()
        self.assertAlmostEqual(self.my_elem.taille_t_plus_dt, 1.25e-03)

    def test_calculer_nouvo_densite(self):
        """ Test de la méthode Element1d.calculer_nouvo_densite() """
        self.my_elem.calculer_nouvo_taille()
        self.my_elem.calculer_nouvo_densite()
        self.assertAlmostEqual(self.my_elem.rho_t_plus_dt, 16258.0)

    def test_calculer_nouvo_pseudo(self):
        """ Test de la méthode Element1d.calculer_nouvo_pseudo() """
        self.my_elem.calculer_nouvo_taille()
        self.my_elem.calculer_nouvo_densite()
        self.my_elem.calculer_nouvo_pseudo(1.0e-6)
        self.assertEqual(self.my_elem.pseudo, 1505370370.3703721)

    def test_calculer_nouvo_dt(self):
        """ Test de la méthode Element1d.calculer_nouvo_dt() """
        self.my_elem.calculer_nouvo_taille()
        self.my_elem.calculer_nouvo_densite()
        self.my_elem.calculer_nouvo_pseudo(1.0e-6)
        self.my_elem.calculer_nouvo_pression()
        self.my_elem.calculer_nouvo_dt()
        self.assertEqual(self.my_elem.delta_t, 7.1888585672553036e-07)

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
