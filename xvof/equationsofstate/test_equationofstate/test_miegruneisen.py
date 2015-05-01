#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module miegrueneisen
"""
import unittest

import xvof.equationsofstate.miegruneisen as mg


class MieGruneisenTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'MieGruneisen'
    """
    def setUp(self):
        """
        Initialisation des tests
        """
        self.equation_of_state = mg.MieGruneisen()
        self.rho = 9.000000000001e+03
        self.e_int = 1.0e+04

    def test_solve_ve(self):
        """ Test de la méthode MieGruneisen.solve_ve() """
        (pression, dpde, vson) = self.equation_of_state.solve_ve(
            1.0 / self.rho, self.e_int)
        self.assertEqual(pression, 1.6111579720792692e+10)
        self.assertEqual(dpde, 1.3441900000000502e+04)
        self.assertEqual(vson, 4.871932359704581e+03)

    def test_solve_vp(self):
        """ Test de la méthode MieGruneisen.solve_vp() """
        dummy = 0.
        with self.assertRaises(NotImplementedError):
            self.equation_of_state.solve_vp(dummy, dummy)

    def test_solve_vt(self):
        """ Test de la méthode MieGruneisen.solve_vt() """
        dummy = 0.
        with self.assertRaises(NotImplementedError):
            self.equation_of_state.solve_vt(dummy, dummy)

if __name__ == '__main__':
    unittest.main()
