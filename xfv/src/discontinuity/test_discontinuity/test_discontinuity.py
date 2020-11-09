#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# pylint: disable=protected-access
"""
Classe de test du module Discontinuity
:TODO : tester l'unicitï¿½ des discontinuitï¿½s identifiï¿½es par label / mask_in/out
"""

import unittest
import os
import numpy as np
from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.data.data_container import DataContainer
from xfv.src.data.enriched_mass_matrix_props import (LumpMenouillardMassMatrixProps,
                                                     LumpSumMassMatrixProps)


class DiscontinuityTest(unittest.TestCase):
    """
        Test case utilisï¿½ pour test les fonctions du module 'Discontinuity'
        """

    def setUp(self):
        """
        Prepare unittests
        """
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_hydro.json")
        self.test_datacontainer = DataContainer(data_file_path)

        self.mask_in = np.array([True, False, False, False])
        self.mask_out = np.array([False, True, False, False])
        self.my_disc = Discontinuity(0, self.mask_in, self.mask_out, 0.2,
                                     LumpMenouillardMassMatrixProps())

    def tearDown(self):
        DataContainer.clear()
        return super().tearDown()

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
        Checks that the position of disc is between 0 and 1
        """
        with self.assertRaises(ValueError):
            Discontinuity(0, self.mask_in, self.mask_out, 2., LumpSumMassMatrixProps())

    def test_detect_hill_disc_node_masks(self):
        """
        Test that a node cannot belong to mask_in and mask_out at the same time
        """
        with self.assertRaises(ValueError):
            Discontinuity(0, np.array([True, False, False, False]),
                          np.array([True, True, False, False]),
                          0.2, LumpSumMassMatrixProps())

    def test_has_mass_matrix_been_computed(self):
        """
        Test of has_mass_matrix_been_computed
        """
        self.my_disc.has_mass_matrix_been_computed()
        np.testing.assert_equal(self.my_disc.mass_matrix_updated, True)


if __name__ == '__main__':
    unittest.main()
