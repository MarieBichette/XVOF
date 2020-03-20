#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
one_dimension_cell module unit tests
"""
import numpy as np
import unittest
import os
from xvof.src.cell.one_dimension_cell import OneDimensionCell
from xvof.src.data.data_container import DataContainer
from xvof.src.cell.test_cell.test_variables import TestVariables


class OneDimensionCellTest(unittest.TestCase):

    def setUp(self):
        data_file_path = os.path.join(os.path.dirname(__file__), "../../../tests/0_UNITTEST/XDATA_epp.xml")
        self.test_datacontainer = DataContainer(data_file_path)

        self.nbr_cells = 4
        self.my_cells = OneDimensionCell(self.nbr_cells)
        self.my_cells.cell_in_target = np.ones(self.nbr_cells, dtype='bool')

    def tearDown(self):
        DataContainer.clear()
        pass

    def test_apply_plastic_corrector_on_deviatoric_stress_tensor(self):
        """
        Test de la méthode apply_plastic_corrector_on_deviatoric_stress_tensor
        """
        mask = np.array([True, True, True, False])
        self.my_cells._deviatoric_stress_new = np.array([[8e+9, -4e+9, -4e+9],
                                                         [5e+8, -2.5e+8, -2.5e+8],
                                                         [1.e+2, -5e+1, -5e+1],
                                                         [4e+9, -2e+9, -2e+9]])
        Y = self.test_datacontainer.material_target.initial_values.yield_stress_init
        self.my_cells.yield_stress.current_value = np.ones(self.nbr_cells) * Y
        self.my_cells.apply_plastic_corrector_on_deviatoric_stress_tensor(mask)
        expected_value =  np.array([[8.000000e+07,  -4.000000e+07,  -4.000000e+07],
                                    [8.000000e+07,  -4.000000e+07,  -4.000000e+07],
                                    [1.e+2, -5e+1, -5e+1],
                                    [4e+9, -2e+9, -2e+9]])
        np.testing.assert_allclose(self.my_cells.deviatoric_stress_new, expected_value)

    def test_compute_equivalent_plastic_strain_rate(self):
        """
        Test de la méthode compute_equivalent_plastic_strain_rate
        """
        mask = np.array([True, True, True, False])
        dt = 1.
        self.my_cells._deviatoric_stress_new = np.array([[8e+9, -4e+9, -4e+9],
                                                             [5e+8, -2.5e+8, -2.5e+8],
                                                             [1.e+2, -5e+1, -5e+1],
                                                             [4e+9, -2e+9, -2e+9]])
        G = self.test_datacontainer.material_target.initial_values.shear_modulus_init
        Y = self.test_datacontainer.material_target.initial_values.yield_stress_init
        self.my_cells.shear_modulus.current_value = np.ones(self.nbr_cells) * G
        self.my_cells.yield_stress.current_value = np.ones(self.nbr_cells) * Y
        self.my_cells.compute_equivalent_plastic_strain_rate(mask, dt)
        np.testing.assert_allclose(self.my_cells.equivalent_plastic_strain_rate,
                                   np.array([0.08301887,  0.00440252,  0.,  0.]), rtol=1.e-5)


if __name__ == "__main__":
    unittest.main()
