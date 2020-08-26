# -*- coding: utf-8 -*-
# pylint: disable=protected-access
"""
Cell module unit tests
"""
import unittest
import unittest.mock as mock
import os
import numpy as np

from xfv.src.cell.one_dimension_enriched_cell_hansbo import OneDimensionHansboEnrichedCell
from xfv.src.cell.one_dimension_cell import OneDimensionCell
from xfv.src.data.data_container import DataContainer
from xfv.src.discontinuity.discontinuity import Discontinuity


class OneDimensionEnrichedHansboCellHydroTest(unittest.TestCase):
    """
    Unittests of the OneDimensionHansboEnrichedCell class
    """

    @classmethod
    def setUpClass(cls):
        """
        Tests setup for class
        """
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_enrichment_hydro.json")
        DataContainer(data_file_path)

    @classmethod
    def tearDownClass(cls):
        DataContainer.clear()
        print("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")

    def setUp(self):
        """
        Unit tests setup
        """
        self.my_cells = OneDimensionHansboEnrichedCell(1)
        self.my_cells._classical = np.array([False])

        self.my_cells.energy.current_value = np.array([1.e+06])
        self.my_cells.pressure.current_value = np.array([1.5e+09])
        self.my_cells.density.current_value = np.array([8000.])
        self.my_cells.pseudo.current_value = np.array([1.e+08])
        self.my_cells.sound_velocity.current_value = np.array([300.])
        self.my_cells._deviatoric_stress_current = np.array([[5., 8., 3.]])

        self.my_cells.energy.new_value = np.array([0.8e+06])
        self.my_cells.pressure.new_value = np.array([1.3e+09])
        self.my_cells.density.new_value = np.array([8020.])
        self.my_cells.pseudo.new_value = np.array([1.e+08])
        self.my_cells.sound_velocity.new_value = np.array([302.])
        self.my_cells._deviatoric_stress_new = np.array([[4., 5., 6.]])
        self.my_cells._deviatoric_strain_rate = np.array([[1., 1., 1.]])

        self.my_cells.enr_density.current_value = np.array([4000.])
        self.my_cells.enr_density.new_value = np.array([4020.])
        self.my_cells.enr_pressure.current_value = np.array([1.1e+09])
        self.my_cells.enr_pressure.new_value = np.array([1.3e+09])
        self.my_cells.enr_energy.current_value = np.array([1.e+06])
        self.my_cells.enr_energy.new_value = np.array([0.8e+06])
        self.my_cells.enr_artificial_viscosity.current_value = np.array([1.e+08])
        self.my_cells.enr_artificial_viscosity.new_value = np.array([1.e+08])
        self.my_cells.enr_sound_velocity.current_value = np.array([300.])
        self.my_cells.enr_sound_velocity.new_value = np.array([302.])
        self.my_cells._enr_deviatoric_stress_current = np.array([[3., 2., 1.], ])
        self.my_cells._enr_deviatoric_stress_new = np.array([[5., 12., 7.], ])
        self.my_cells._enr_deviatoric_stress_new = np.array([[5., 12., 7.], ])
        self.my_cells._enr_deviatoric_strain_rate = np.array([[4., 3., 8.], ])
        self.my_cells.enr_yield_stress.current_value = np.array([10.])
        self.my_cells._enr_equivalent_plastic_strain_rate = np.array([0.])
        self.my_cells._enr_stress = np.array([[0., 0., 0.]])
        self.my_cells.left_part_size.current_value = np.array([0.2])
        self.my_cells.right_part_size.current_value = np.array([0.3])
        self.my_cells.left_part_size.new_value = np.array([0.4])
        self.my_cells.right_part_size.new_value = np.array([0.6])

        # configuration d'un mock 'discontinuity'
        config = {'mask_in_nodes': np.array([True, False]),
                  'mask_out_nodes': np.array([False, True]),
                  'position_in_ruptured_element': (
                      DataContainer().material_target.failure_model.failure_treatment_value),
                  'mask_ruptured_cell': np.array([True]),
                  'ruptured_cell_id': np.array([0]),
                  'plastic_cells': np.array([False]),
                  'enr_velocity_new': np.array([[1., ], [3., ]])
                }
        self.__patcher = mock.patch('xfv.src.discontinuity.discontinuity.Discontinuity',
                             spec=Discontinuity, **config)
        self.mock_discontinuity = self.__patcher.start()

    def tearDown(self):
        self.__patcher.stop()
        return super().tearDown()

    def test_classical(self):
        """
        Test of the property classical
        """
        np.testing.assert_array_equal(self.my_cells.classical, np.zeros([1], dtype="bool"))

    def test_enriched(self):
        """
        Test of the property enriched
        """
        np.testing.assert_array_equal(self.my_cells.enriched, np.ones([1], dtype="bool"))

    @mock.patch.object(OneDimensionCell, "apply_equation_of_state",
                       spec=classmethod, new_callable=mock.MagicMock)
    @mock.patch.object(OneDimensionCell, "add_elastic_energy_method",
                       spec=classmethod, new_callable=mock.MagicMock)
    def test_compute_enriched_elements_new_pressure_without_elasticity(
            self, mock_elasto, mock_eos):  #pylint: disable=unused-argument
        """
        Test of the compute_enriched_elements_new_pressure method (Hansbo case)
        """
        # Configuration des mocks
        mock_eos.side_effect = [[self.my_cells.enr_energy.new_value,
                                 self.my_cells.enr_pressure.new_value,
                                 self.my_cells.enr_sound_velocity.new_value]]

        self.my_cells.compute_enriched_elements_new_pressure(1.)

        mock_eos.assert_any_call(
            self.my_cells, self.my_cells._target_eos,
            self.my_cells.enr_density.current_value,
            self.my_cells.enr_density.new_value,
            self.my_cells.enr_pressure.current_value,
            self.my_cells.enr_pressure.new_value,
            self.my_cells.enr_energy.current_value,
            self.my_cells.enr_energy.new_value,
            self.my_cells.enr_artificial_viscosity.current_value,
            self.my_cells.enr_sound_velocity.new_value)


if __name__ == "__main__":
    unittest.main()
