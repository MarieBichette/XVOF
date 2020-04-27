#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# pylint: disable=protected-access
"""
one_dimension_cell module unit tests
"""
import unittest
import os
import numpy as np
from xfv.src.cell.one_dimension_cell import OneDimensionCell
from xfv.src.data.data_container import DataContainer


class OneDimensionCellEPPTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Tests setup for class
        """
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_epp.json")
        DataContainer(data_file_path)

    @classmethod
    def tearDownClass(cls):
        DataContainer.clear()
        print("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")

    def setUp(self):
        self.nbr_cells = 4
        self.my_cells = OneDimensionCell(self.nbr_cells)
        self.my_cells.cell_in_target = np.ones(self.nbr_cells, dtype='bool')
        self.test_data = DataContainer()  # pylint: disable=no-value-for-parameter

    def tearDown(self):
        pass

    def test_apply_plastic_corrector_on_deviatoric_stress_tensor(self):
        """
        Test de la mï¿½thode apply_plastic_corrector_on_deviatoric_stress_tensor
        """
        mask = np.array([True, True, False, False])
        self.my_cells._deviatoric_stress_new = np.array([[8e+9, -4e+9, -4e+9],
                                                         [5e+8, -2.5e+8, -2.5e+8],
                                                         [1.e+2, -5e+1, -5e+1],
                                                         [4e+9, -2e+9, -2e+9]])
        yield_stress = self.test_data.material_target.initial_values.yield_stress_init
        self.my_cells.yield_stress.current_value = np.ones(self.nbr_cells) * yield_stress
        self.my_cells.apply_plastic_corrector_on_deviatoric_stress_tensor(mask)
        expected_value = np.array([[8.000000e+07,  -4.000000e+07,  -4.000000e+07],
                                   [8.000000e+07,  -4.000000e+07,  -4.000000e+07],
                                   [1.e+2, -5e+1, -5e+1],
                                   [4e+9, -2e+9, -2e+9]])
        np.testing.assert_allclose(self.my_cells.deviatoric_stress_new, expected_value)

    def test_compute_equivalent_plastic_strain_rate(self):
        """
        Test de la mï¿½thode compute_equivalent_plastic_strain_rate
        """
        mask = np.array([True, True, False, False])
        delta_t = 1.
        self.my_cells._deviatoric_stress_new = np.array([[8e+9, -4e+9, -4e+9],
                                                         [5e+8, -2.5e+8, -2.5e+8],
                                                         [1.e+2, -5e+1, -5e+1],
                                                         [4e+9, -2e+9, -2e+9]])
        shear_modulus = self.test_data.material_target.initial_values.shear_modulus_init
        yield_stress = self.test_data.material_target.initial_values.yield_stress_init
        self.my_cells.shear_modulus.current_value = np.ones(self.nbr_cells) * shear_modulus
        self.my_cells.yield_stress.current_value = np.ones(self.nbr_cells) * yield_stress
        self.my_cells.compute_equivalent_plastic_strain_rate(mask, delta_t)
        np.testing.assert_allclose(self.my_cells.equivalent_plastic_strain_rate,
                                   np.array([0.08301887,  0.00440252,  0.,  0.]), rtol=1.e-5)


if __name__ == "__main__":
    unittest.main()
