#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module miegrueneisen
"""
import xvof.equationsofstate.miegruneisen as mg
import unittest


class MieGruneisenTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'MieGruneisen'
    """
    def setUp(self):
        """
        Initialisation des tests
        """
        self.equation_of_state = mg.MieGruneisen()
        self.rho = 9.000001000003e+03
        self.e_int = 2.0e+03

    def test_solve_ve(self):
        """
        Test de la fonction MieGruneisen.solve_ve
        """
        (pression, dpde, vson) = self.equation_of_state.solve_ve(
            1.0 / self.rho, self.e_int)
        self.assertEqual(pression, 16004065578.856352)
        self.assertEqual(dpde, 13441.9005000015)
        self.assertEqual(vson, 4869.690880710129)

    def test_solve_vp(self):
        """
        Test de la fonction MieGruneisen.solve_vp
        """
        dummy = 0.
        with self.assertRaises(NotImplementedError):
            self.equation_of_state.solve_vp(dummy, dummy)

    def test_solve_vt(self):
        """
        Test de la fonction MieGruneisen.solve_vt
        """
        dummy = 0.
        with self.assertRaises(NotImplementedError):
            self.equation_of_state.solve_vt(dummy, dummy)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############ PROGRAMME PRINCIPAL ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
if __name__ == '__main__':
    unittest.main()