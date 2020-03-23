#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Cell module unit tests
"""
import numpy as np
import unittest
import mock
import os
from xfv.src.mesh.topology1d import Topology1D
from xfv.src.cell.one_dimension_enriched_cell_Hansbo import OneDimensionHansboEnrichedCell
from xfv.src.cell.one_dimension_enriched_cell import OneDimensionEnrichedCell
from xfv.src.cell.one_dimension_cell import OneDimensionCell
from xfv.src.data.data_container import DataContainer
from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.rheology.constantshearmodulus import ConstantShearModulus


class OneDimensionEnrichedHansboCellHydroTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        data_file_path = os.path.join(os.path.dirname(__file__), "../../../tests/0_UNITTEST/XDATA_enrichment_hydro.xml")
        DataContainer(data_file_path)

    @classmethod
    def tearDownClass(cls):
        DataContainer.clear()
        print("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
        pass

    def setUp(self):
        """
        Pr�paration des tests
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

        # configuration d'un mock 'discontinuity'
        config = {'mask_in_nodes': np.array([True, False]),
                  'mask_out_nodes': np.array([False, True]),
                  'position_in_ruptured_element': DataContainer().material_target.failure_model.failure_treatment_value,
                  'mask_ruptured_cell': np.array([True]),
                  'ruptured_cell_id': np.array([0]),
                  'plastic_cells': False,
                  'additional_dof_density.current_value':np.array([4000.]),
                  'additional_dof_density.new_value': np.array([4020.]),
                  'additional_dof_pressure.current_value': np.array([1.1e+09]),
                  'additional_dof_pressure.new_value': np.array([1.3e+09]),
                  'additional_dof_energy.current_value': np.array([1.e+06]),
                  'additional_dof_energy.new_value': np.array([0.8e+06]),
                  'additional_dof_artificial_viscosity.current_value': np.array([1.e+08]),
                  'additional_dof_artificial_viscosity.new_value': np.array([1.e+08]),
                  'additional_dof_sound_velocity.current_value': np.array([300.]),
                  'additional_dof_sound_velocity.new_value': np.array([302.]),
                  'additional_dof_deviatoric_stress_current': np.array([[3., 2., 1.],]),
                  'additional_dof_deviatoric_stress_new': np.array([[5., 12., 7.],]),
                  '_additional_dof_deviatoric_stress_new': np.array([[5., 12., 7.], ]),
                  'additional_dof_deviatoric_strain_rate': np.array([[4., 3., 8.],]),
                  'additional_dof_yield_stress.current_value': np.array([10.]),
                  '_additional_dof_equivalent_plastic_strain_rate': np.array([0.]),
                  'additional_dof_velocity_new': np.array([[1., ], [3., ]]),
                  'additional_dof_stress': np.array([[0., 0., 0.]]),
                  'left_part_size.current_value': np.array([0.2]),
                  'right_part_size.current_value': np.array([0.3]),
                  'left_part_size.new_value': np.array([0.4]),
                  'right_part_size.new_value': np.array([0.6]),
                  }
        patcher = mock.patch('xfv.src.discontinuity.discontinuity.Discontinuity', spec=Discontinuity, **config)
        self.mock_discontinuity = patcher.start()

    def tearDown(self):
        pass

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    @mock.patch.object(OneDimensionCell, "apply_equation_of_state", spec=classmethod, new_callable=mock.MagicMock)
    @mock.patch.object(OneDimensionCell, "add_elastic_energy_method", spec=classmethod, new_callable=mock.MagicMock)
    def test_compute_enriched_elements_new_pressure_without_elasticity(self, mock_elasto, mock_eos, mock_disc_list):
        """
        Test de la m�thode compute_enriched_elements_new_pressure pour Hansbo
        """
        # Configuration des mocks
        Discontinuity.discontinuity_list.return_value = [self.mock_discontinuity]
        mock_eos.side_effect = [[self.my_cells.energy.new_value,
                                 self.my_cells.pressure.new_value,
                                 self.my_cells.sound_velocity.new_value],
                                [self.mock_discontinuity.additional_dof_energy.new_value,
                                 self.mock_discontinuity.additional_dof_pressure.new_value,
                                 self.mock_discontinuity.additional_dof_sound_velocity.new_value]]
        # mock_eos.return_value =

        self.my_cells.compute_enriched_elements_new_pressure(1.)

        mock_elasto.assert_not_called()

        mock_eos.assert_any_call(self.my_cells, DataContainer().material_target.constitutive_model.eos,
                                self.my_cells.density.current_value, self.my_cells.density.new_value,
                                self.my_cells.pressure.current_value, self.my_cells.pressure.new_value,
                                self.my_cells.energy.current_value, self.my_cells.energy.new_value,
                                self.my_cells.pseudo.current_value, self.my_cells.sound_velocity.new_value)

        mock_eos.assert_any_call(self.my_cells, DataContainer().material_target.constitutive_model.eos,
                                self.mock_discontinuity.additional_dof_density.current_value,
                                self.mock_discontinuity.additional_dof_density.new_value,
                                self.mock_discontinuity.additional_dof_pressure.current_value,
                                self.mock_discontinuity.additional_dof_pressure.new_value,
                                self.mock_discontinuity.additional_dof_energy.current_value,
                                self.mock_discontinuity.additional_dof_energy.new_value,
                                self.mock_discontinuity.additional_dof_artificial_viscosity.current_value,
                                self.mock_discontinuity.additional_dof_sound_velocity.new_value)


if __name__ == "__main__":
    unittest.main()
