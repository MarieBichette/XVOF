#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module Discontinuity
:TODO : tester l'unicité des discontinuités identifiées par label / mask_in/out
"""

import unittest
import numpy as np
from xfv.src.mass_matrix.mass_matrix_utilities import inverse_masse, lump_matrix


class MatrixUtliitiesTest(unittest.TestCase):
    """
        Test case utilisé pour test les fonctions du module 'Discontinuity'
        """
    def setUp(self):
        """
        Préparation des tests unitaires
        """
        self._test_matrix = np.array([[1., 2.], [1., 4.]])
        self._test_vecteur = np.array([1., 2.])

    def test_inverse_masse(self):
        """
        Test de la méthode inverse_masse
        """
        # cas matrice = vecteur
        inverse = inverse_masse(self._test_matrix)
        inverse_solution = np.array([[2, -1], [-0.5, 0.5]])
        np.testing.assert_array_equal(inverse, inverse_solution)
        # cas matrice = vraie matrice
        inverse = inverse_masse(self._test_vecteur)
        inverse_solution = np.array([[1.], [0.5]])
        np.testing.assert_array_equal(inverse, inverse_solution)

    def test_lump_matrix(self):
        """
        Test de la méthode lump_matrix
        """
        lump = lump_matrix(self._test_matrix)
        solution = np.array([[3., 0.], [0., 5.]])
        np.testing.assert_array_equal(lump, solution)


if __name__ == '__main__':
    unittest.main()
