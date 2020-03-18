#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
one_dimension_cell module unit tests
"""
import numpy as np
import unittest
import mock
import os
from xvof.src.cell.one_dimension_cell import OneDimensionCell
from xvof.src.mesh.topology1d import Topology1D
from xvof.src.data.data_container import DataContainer
from xvof.src.cell.test_cell.test_variables import TestVariables


class OneDimensionCellTest(unittest.TestCase):

    def setUp(self):
        data_file_path = os.path.join(os.path.dirname(__file__), "../../../tests/0_UNITTEST/XDATA_hydro.xml")
        self.test_datacontainer = DataContainer(data_file_path)
        self.nbr_cells = 4
        self.my_cells = OneDimensionCell(self.nbr_cells)
        self.my_cells.cell_in_target = np.ones(self.nbr_cells, dtype='bool')

    def tearDown(self):
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        pass

    def test_compute_size(self):
        """
        Test of compute_size method
        """
        topo = mock.MagicMock(Topology1D)
        type(topo).nodes_belonging_to_cell = mock.PropertyMock(return_value=np.array([[0, 1], [1, 2], [2, 3], [3, 4]]))
        node_coord = np.array([[-0.35, ], [0.15, ], [0.2, ], [0.25, ], [0.5, ]])
        self.my_cells.compute_size(topo, node_coord)
        np.testing.assert_allclose(self.my_cells.size_t, np.array([0.5, 0.05, 0.05, 0.25]))

    def test_compute_new_size(self):
        """
        Test of compute_new_size method
        """
        topo = mock.MagicMock(Topology1D)
        type(topo).nodes_belonging_to_cell = mock.PropertyMock(return_value=np.array([[0, 1], [1, 2], [2, 3],[3, 4]]))
        mask = np.array([True, True, False, False])
        self.my_cells._size_t = np.array([0.5, 0.05, 0.05, 0.25])
        node_new_coord = np.array([[-0.25, ], [0.1, ], [0.2, ], [0.45, ], [0.85, ]])
        self.my_cells.compute_new_size(topo, node_new_coord, mask)
        np.testing.assert_allclose(self.my_cells.size_t_plus_dt[mask], np.array([0.35, 0.1, 0.25, 0.4]))

    @mock.patch.object(geometrical_props, "section", new_callable=mock.PropertyMock, return_value=0.0003141592653589793)
    def test_compute_mass(self, mock_geom_props):
        """
        Test of compute_mass method
        """
        self.my_cells.density.current_value = np.array([8129., 8129., 8129., 8129.])
        self.my_cells._size_t = np.array([0.6, 0.1, 0.15, 0.6])
        self.test_cell.compute_mass()
        np.testing.assert_allclose(self.test_cell.mass, [1.5322804 , 0.25538007, 0.3830701, 1.5322804])

    def test_compute_new_density(self):
        """
        Test of compute_new_density method
        """
        self.my_cells.density.current_value = np.ones(self.nb_cells) * 8930.
        self.my_cells._size_t = np.array([0.6, 0.1, 0.15])
        self.my_cells._size_t_plus_dt = np.array([0.8, 0.1, 0.27])
        mask = np.array([True, True, True, False])
        self.my_cells.compute_new_density(mask)
        np.testing.assert_allclose(self.my_cells.density.new_value, np.array([6096.75, 8129., 4516.11111111, 8930.]))

    def test_compute_new_pressure_without_elasticity(self):
        """
        Test de la méthode compute_new_pressure
        """
        mask = np.array([True, True, True, False])
        dt = 1.

        self.my_cells.density.current_value = np.copy(self.test_variables.density_old)
        self.my_cells.pressure.current_value = np.copy(self.test_variables.pressure_old)
        self.my_cells.energy.current_value = np.copy(self.test_variables.energy_old)

        self.my_cells.density.new_value = np.copy(self.test_variables.density_new)
        self.my_cells.sound_velocity.new_value = np.zeros([self.nbr_cells])
        self.my_cells.pressure.new_value = np.zeros([self.nbr_cells])
        self.my_cells.energy.new_value = np.zeros([self.nbr_cells])
        self.my_cells.pseudo.new_value = np.zeros([self.nbr_cells])

        self.my_cells.compute_new_pressure(mask, dt)
        mock_add_elasticity.assert_not_called()

        np.testing.assert_allclose(self.my_cells.energy.new_value[mask], self.test_variables.energy_new[mask])
        np.testing.assert_allclose(self.my_cells.energy.new_value[~mask], self.test_variables.energy_old[~mask])
        np.testing.assert_allclose(self.my_cells.pressure.new_value[mask], self.test_variables.pressure_new[mask])
        np.testing.assert_allclose(self.my_cells.pressure.new_value[~mask], self.test_variables.pressure_old[~mask])

    def test_compute_new_pseudo(self):
        """
        Test of compute_new_pseudo
        """
        mask = [True, True, False, False]
        self.my_cells.density.current_value = np.array([8000., 8500., 9500., 8000.])
        self.my_cells.density.new_value = np.array([8120., 8440., 9620., 8120.])
        self.my_cells.sound_velocity.current_value = np.array([300., 300., 300., 300.])
        self.my_cells.pseudo.new_value = np.array([1.e+9, 2.e+9, 3.e+9, 1.e+9])
        self.my_cells._size_t_plus_dt = np.array([1., 1., 2., 1.])
        assert self.test_datacontainer.numeric.a_pseudo == 2.
        assert self.test_datacontainer.numeric.b_pseudo == 0.2
        self.my_cells.pseudo.new_value[mask] = self.my_cells.compute_new_pseudo(2., mask)
        np.testing.assert_allclose(self.my_cells.pseudo.new_value[mask],
                                   np.array([0.446760955365777, 0.2125398482595006, 3.e+9, 1.e+9]))

    def test_compute_new_time_step(self):
        """
        Test de la méthode compute_enriched_elements_new_time_step
        """
        self.my_cells.density.current_value = np.array([8000., 8500., 9500., 8000.])
        self.my_cells.density.new_value = np.array([8120., 8440., 9000., 8120.])
        self.my_cells.sound_velocity.new_value = np.array([100., 100., 100., 100.])
        self.my_cells.pseudo.current_value = np.array([0., 0., 0., 0.])
        self.my_cells.pseudo.new_value = np.array([1., 1., 1., 1.])
        self.my_cells._size_t_plus_dt = np.array([1., 1., 2., 1.])
        assert self.test_datacontainer.numeric.cfl == 0.34
        assert self.test_datacontainer.numeric.cfl_pseudo == 0.
        self.my_cells._dt = np.array([1., 2., 5., 6.])
        mask = np.array([False, True, True, False])
        self.my_cells.compute_new_time_step(mask)
        np.testing.assert_allclose(self.my_cells._dt, np.array([1., 3.400000e-03, 6.800000e-03, 6.]))

    def test_compute_complete_stress_tensor_hydro(self):
        """
        Test de la méthode compute_complete_stress_tensor
        """
        mask = self.mask
        self.my_cells.energy.current_value = np.copy(self.test_variables.energy_old)
        self.my_cells._deviatoric_strain_rate = np.copy(self.test_variables.deviatoric_strain_rate)
        self.my_cells._deviatoric_stress_new = np.copy(self.test_variables.deviatoric_stress_new)
        self.my_cells._deviatoric_stress_current = np.copy(self.test_variables.deviatoric_stress_current)
        self.my_cells.pressure.new_value = np.copy(self.test_variables.pressure_new)
        self.my_cells.pseudo.new_value = np.copy(self.test_variables.pseudo_new)
        self.my_cells.compute_complete_stress_tensor(mask)
        np.testing.assert_allclose(self.my_cells.stress[mask], self.test_variables.stress_new[mask])

    def test_impose_pressure(self):
        """
        Test de impose_pressure
        """
        self.my_cells.pressure.new_value = np.copy(self.test_variables.pressure_new)
        self.my_cells.impose_pressure(0, -1.)
        np.testing.assert_array_equal(self.my_cells.pressure.new_value[1:], self.test_variables.pressure_new[1:])
        np.testing.assert_array_equal(self.my_cells.pressure.new_value[0], np.array([-1.]))

if __name__ == "__main__":
    unittest.main()
