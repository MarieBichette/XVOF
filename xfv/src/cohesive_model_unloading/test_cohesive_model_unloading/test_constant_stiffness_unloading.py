# -*- coding: utf-8 -*-
"""
Class to test the module ConstantStiffnessUnloadingTest
"""
import unittest
import os
import numpy as np
from xfv.src.cohesive_model_unloading.constant_stiffness_unloading import ConstantStiffnessUnloading
from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.data.data_container import DataContainer
from xfv.src.data.enriched_mass_matrix_props import LumpMenouillardMassMatrixProps


class ConstantStiffnessUnloadingTest(unittest.TestCase):
    """
    Test case for ConstantStiffnessUnloadingTest
    """
    def setUp(self):
        """
        Tests initialisation
        """
        # DataContainer creation (just to be able to build the discontinuity)
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_hydro.json")
        self.test_datacontainer = DataContainer(data_file_path)
        # Discontinuity creation
        self.disc = Discontinuity(np.array([True, False]), np.array([False, True]), 0.5,
                                  LumpMenouillardMassMatrixProps())
        # Creation of the tested service
        self.test_unloading_model = ConstantStiffnessUnloading(10.)

    def test_compute_unloading_reloading_condition(self):
        """
        Test of the method compute_unloading_reloading_condition du module
        ConstantStiffnessUnloading
        """
        self.disc.history_min_cohesive_force = 40.
        self.disc.history_max_opening = 2.
        result = self.test_unloading_model.compute_unloading_reloading_condition(self.disc, 0.5)
        expected = 25.
        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
