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
        data_file_path = os.path.realpath(os.path.join(os.getcwd(), "../tests/0_UNITTEST/XDATA_epp.xml"))
        self.test_datacontainer = DataContainer(data_file_path)
        self.nbr_cells = 4
        self.my_cells = OneDimensionCell(self.nbr_cells)
        self.my_cells.cell_in_target = np.ones(self.nbr_cells, dtype='bool')

        self.test_variables = TestVariables(4, 5)
        self.test_variables.define_epp_variables()

        self.mask = np.array([True, True, False, False])

    def tearDown(self):
        pass

    def test_apply_plastic_corrector_on_deviatoric_stress_tensor(self):
        """
        Test de la méthode apply_plastic_corrector_on_deviatoric_stress_tensor
        """
        print __name__ + " : Test apply correction on deviatoric stress tensor"
        mask = self.mask
        self.my_cells._deviatoric_stress_new = np.copy(self.test_variables.deviatoric_stress_new)
        self.my_cells.yield_stress.current_value = np.copy(self.test_variables.yield_stress_old)
        self.my_cells.apply_plastic_corrector_on_deviatoric_stress_tensor(mask)
        np.testing.assert_allclose(self.my_cells.deviatoric_stress_new,
                                   self.test_variables.deviatoric_stress_new_after_plastic_correction)
        print "__[OK]"

    def test_compute_equivalent_plastic_strain_rate(self):
        """
        Test de la méthode compute_equivalent_plastic_strain_rate
        """
        print __name__ + " : Test compute equivalent plastic strain rate"
        mask = self.mask
        dt = self.test_variables.dt
        self.my_cells._deviatoric_stress_new =np.copy(self.test_variables.deviatoric_stress_new)
        self.my_cells.shear_modulus.current_value = np.copy(self.test_variables.shear_modulus_old)
        self.my_cells.yield_stress.current_value = np.copy(self.test_variables.yield_stress_old)
        self.my_cells.compute_equivalent_plastic_strain_rate(mask, dt)
        np.testing.assert_allclose(self.my_cells.equivalent_plastic_strain_rate[mask],
                                   self.test_variables.equivalent_plastic_strain_rate[mask])
        print "__[OK]"

if __name__ == "__main__":
    unittest.main()
