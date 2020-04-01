#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module Discontinuity
:TODO : tester l'unicit� des discontinuit�s identifi�es par label / mask_in/out
"""

import unittest
import numpy as np
import os
from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.data.data_container import DataContainer


class DiscontinuityTest(unittest.TestCase):
    """
        Test case utilis� pour test les fonctions du module 'Discontinuity'
        """

    def setUp(self):
        """
        Pr�paration des tests unitaires
        """
        data_file_path = os.path.join(os.path.dirname(__file__), "../../../tests/0_UNITTEST/XDATA_hydro.json")
        self.test_datacontainer = DataContainer(data_file_path)

        self.mask_in = np.array([True, False, False, False])
        self.mask_out = np.array([False, True, False, False])
        self.my_disc = Discontinuity(self.mask_in, self.mask_out, 0.2)

    def tearDown(self):
        DataContainer.clear()
        return super(DiscontinuityTest, self).tearDown()


    # def test_label(self):
    #     """
    #     Teste l'accesseur label de la classe Discontinuity
    #     """
    #     np.testing.assert_equal(self.my_disc.label, 1.)

    def test_mask_in_nodes(self):
        """
        Teste l'accesseur mask_in_nodes de la classe Discontinuity
        """
        np.testing.assert_equal(self.my_disc.mask_in_nodes, np.array([True, False, False, False]))

    def test_mask_out_nodes(self):
        """
        Teste l'accesseur mask_out_nodes de la classe Discontinuity
        """
        np.testing.assert_equal(self.my_disc.mask_out_nodes, np.array([False, True, False, False]))


    def test_detect_hill_disc_position(self):
        """
        V�rifie que la position de la disc est comprise entre 0 et 1
        """
        with self.assertRaises(ValueError):
            Discontinuity(self.mask_in, self.mask_out, 2.)

    def test_detect_hill_disc_node_masks(self):
        """
        Test qu'un noeud ne peut pas appartenir au mask in et out en m�me temps
        """
        with self.assertRaises(ValueError):
            Discontinuity(np.array([True, False, False, False]), np.array([True, True, False, False]), 0.2)


    # def test_mass_matrix_updated_true(self):
    #     """Test de la m�thode mass matrix updated """
    #     self.my_disc.__mass_matrix_updated = True # marche pas, ni comme �a, ni avec mock
    #     np.testing.assert_equal(self.my_disc.mass_matrix_updated, True)
    #     self.my_disc.__mass_matrix_updated = False
    #     np.testing.assert_equal(self.my_disc.mass_matrix_updated, False)

    def test_hasMassMatrixBeenComputed(self):
        """Teste de la m�thode has_mass_matrix_been_computed pour le module Discontinuity"""
        # utilise la propri�t� mass_matrix_updated test�e plus haut
        self.my_disc.has_mass_matrix_been_computed()
        np.testing.assert_equal(self.my_disc.mass_matrix_updated, True)

if __name__ == '__main__':
    unittest.main()
