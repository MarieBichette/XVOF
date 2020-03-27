#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Classe de test du module cohesive unloading
"""
import unittest
import numpy as np
import os
from xfv.src.cohesive_model_unloading.loss_of_stiffness_unloading import LossOfStiffnessUnloading
from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.data.data_container import DataContainer


class LossOfStiffnessUnloadingTest(unittest.TestCase):
    """
    Test case for ConstantStiffnessUnloadingTest
    """
    def setUp(self):
        """
        Initialisation des tests
        """
        # DataContainer creation (just to be able to build the discontinuity)
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_hydro.xml")
        self.test_datacontainer = DataContainer(data_file_path)
        # Discontinuity creation
        self.disc = Discontinuity(np.array([True, False]), np.array([False, True]), 0.5)
        # Creation of the tested service
        self.test_unloading_model = LossOfStiffnessUnloading()

    def tearDown(self):
        """
        Operations to be done after completing all the tests in the class
        """
        pass

    def test_compute_unloading_reloading_condition(self):
        """
        Test of the method compute_unloading_reloading_condition du module
        ConstantStiffnessUnloading
        """
        self.disc.history_min_cohesive_force = 40.
        self.disc.history_max_opening = 2.
        result = self.test_unloading_model.compute_unloading_reloading_condition(self.disc, 0.5)
        expected = 10.
        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
