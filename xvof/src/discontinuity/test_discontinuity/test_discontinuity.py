#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module Discontinuity
:TODO : tester l'unicité des discontinuités identifiées par label / mask_in/out
"""

import unittest
import numpy as np
from xvof.src.discontinuity.discontinuity import Discontinuity


class DiscontinuityTest(unittest.TestCase):
    """
        Test case utilisé pour test les fonctions du module 'Discontinuity'
        """

    def setUp(self):
        """
        Préparation des tests unitaires
        """
        self.mask_in = np.array([True, False, False, False])
        self.mask_out = np.array([False, True, False, False])
        self.my_disc = Discontinuity(self.mask_in, self.mask_out)

    # def test_label(self):
    #     """
    #     Teste l'accesseur label de la classe Discontinuity
    #     """
    #     np.testing.assert_equal(self.my_disc.label, 1.)
    #
    # def test_mask_in_nodes(self):
    #     """
    #     Teste l'accesseur mask_in_nodes de la classe Discontinuity
    #     """
    #     np.testing.assert_equal(self.my_disc.mask_in_nodes, np.array([True, False, False, False]))
    #
    # def test_mask_out_nodes(self):
    #     """
    #     Teste l'accesseur mask_out_nodes de la classe Discontinuity
    #     """
    #     np.testing.assert_equal(self.my_disc.mask_out_nodes, np.array([False, True, False, False]))
    #
    # def test_mass_matrix_updated_true(self):
    #     """Test de la méthode mass matrix updated """
    #     self.my_disc.__mass_matrix_updated = True # marche pas, ni comme ça, ni avec mock
    #     np.testing.assert_equal(self.my_disc.mass_matrix_updated, True)
    #     self.my_disc.__mass_matrix_updated = False
    #     np.testing.assert_equal(self.my_disc.mass_matrix_updated, False)

    def test_hasMassMatrixBeenComputed(self):
        """Teste de la méthode has_mass_matrix_been_computed pour le module Discontinuity"""
        # utilise la propriété mass_matrix_updated testée plus haut
        self.my_disc.has_mass_matrix_been_computed()
        np.testing.assert_equal(self.my_disc.mass_matrix_updated, True)

if __name__ == '__main__':
    unittest.main()
