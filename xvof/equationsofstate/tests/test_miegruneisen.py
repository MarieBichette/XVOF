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

    def testsolveVolumeEnergy(self):
        """ Test de la méthode MieGruneisen.solveVolumeEnergy() """
        (pression, dpde, vson) = self.equation_of_state.solveVolumeEnergy(
            1.0 / self.rho, self.e_int)
        self.assertEqual(pression, 1.6111579720792692e+10)
        self.assertEqual(dpde, 1.3441900000000502e+04)
        self.assertEqual(vson, 4.871932359704581e+03)

    def testsolveVolumePressure(self):
        """ Test de la méthode MieGruneisen.solveVolumePressure() """
        dummy = 0.
        with self.assertRaises(NotImplementedError):
            self.equation_of_state.solveVolumePressure(dummy, dummy)

    def testsolveVolumeTemperature(self):
        """ Test de la méthode MieGruneisen.solveVolumeTemperature() """
        dummy = 0.
        with self.assertRaises(NotImplementedError):
            self.equation_of_state.solveVolumeTemperature(dummy, dummy)

if __name__ == '__main__':
    unittest.main()
