#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
one_dimension_cell module unit tests
"""
import numpy as np
import unittest
import mock
from xvof.src.cell.one_dimension_cell import OneDimensionCell
from xvof.src.mesh.topology1d import Topology1D
from xvof.src.data.data_container import DataContainer
from xvof.src.cell.test_cell.test_variables import TestVariables


class OneDimensionCellTest(unittest.TestCase):

    def setUp(self):
        data_file_path = "//home/marie/PycharmProjects/XVOF/xvof.src/0_UNITTEST/XDATA_hydro.xml"
        self.test_datacontainer = DataContainer(data_file_path)
        self.nbr_cells = 4
        self.my_cells = OneDimensionCell(self.nbr_cells)
        self.my_cells.cell_in_target = np.ones(self.nbr_cells, dtype='bool')

        self.test_variables = TestVariables(4, 5)
        self.test_variables.define_hydro_variables()

        self.mask = np.array([True, True, False, False])

    def tearDown(self):
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        pass

    def test_compute_size(self):
        """
        Test of compute_size method
        """
        print "Test " + __name__
        topo = mock.MagicMock(Topology1D)
        type(topo).nodes_belonging_to_cell = mock.PropertyMock(return_value=np.array([[0, 1], [1, 2], [2, 3], [3, 4]]))
        self.my_cells.compute_size(topo, self.test_variables.node_coord_old)
        np.testing.assert_allclose(self.my_cells.size_t, self.test_variables.cell_size_old)
        print "__[OK]"

    def test_compute_new_size(self):
        """
        Test of compute_new_size method
        """
        print "Test " + __name__
        topo = mock.MagicMock(Topology1D)
        type(topo).nodes_belonging_to_cell = mock.PropertyMock(return_value=np.array([[0, 1], [1, 2], [2, 3],[3, 4]]))
        mask = self.mask
        self.my_cells._size_t = np.copy(self.test_variables.cell_size_old)
        self.my_cells.compute_new_size(topo, self.test_variables.node_coord_new, mask)
        np.testing.assert_allclose(self.my_cells.size_t_plus_dt[mask], self.test_variables.cell_size_new[mask])
        np.testing.assert_allclose(self.my_cells.size_t_plus_dt[~mask], np.zeros([self.test_variables.nb_cells])[~mask])
        print "__[OK]"

    def test_compute_mass(self):
        """
        Test of compute_mass method
        """
        print "Test " + __name__
        self.my_cells.density.current_value = np.copy(self.test_variables.density_old)
        self.my_cells._size_t = np.copy(self.test_variables.cell_size_old)
        self.my_cells.compute_mass()
        np.testing.assert_allclose(self.my_cells.mass, self.test_variables.cell_mass)
        print "__[OK]"

    def test_compute_new_density(self):
        """
        Test of compute_new_density method
        """
        print "Test" + __name__
        self.my_cells.density.current_value = np.copy(self.test_variables.density_init)
        self.my_cells._size_t = np.copy(self.test_variables.cell_size_init)
        self.my_cells._size_t_plus_dt = np.copy(self.test_variables.cell_size_old)
        mask = self.mask
        self.my_cells.compute_new_density(mask)
        np.testing.assert_allclose(self.my_cells.density.new_value[mask], self.test_variables.density_old[mask])
        np.testing.assert_allclose(self.my_cells.density.new_value[~mask], self.test_variables.density_init[~mask])
        print "__[OK]"

    @mock.patch.object(OneDimensionCell, "add_elastic_energy_method", spec=classmethod, new_callable=mock.MagicMock)
    def test_compute_new_pressure_without_elasticity(self, mock_add_elasticity):
        """
        Test de la méthode compute_new_pressure
        """
        print "Test " + __name__
        # L'option élasticité est désactivée dans le jeu de donnée
        type(DataContainer().material_target.constitutive_model).elasticity_model = mock.PropertyMock(return_value=None)
        mask = self.mask
        dt = self.test_variables.dt

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
        print "__[OK]"

    def test_compute_new_pseudo(self):
        """
        Test of compute_new_pseudo
        """
        print "Test " + __name__
        mask = self.mask
        dt = self.test_variables.dt
        self.my_cells.pseudo.current_value = np.copy(self.test_variables.pseudo_old)
        self.my_cells.compute_new_pseudo(dt, mask)
        np.testing.assert_allclose(self.my_cells.pseudo.new_value[mask], self.test_variables.pseudo_new[mask])
        np.testing.assert_allclose(self.my_cells.pseudo.new_value[~mask],
                                   np.zeros([self.test_variables.nb_cells])[~mask])
        print "__[OK]"

    @mock.patch.object(OneDimensionCell, "compute_time_step", spec=classmethod, new_callable=mock.MagicMock)
    def test_compute_new_time_step(self, mock_compute_dt):
        """
        Test de la méthode compute_enriched_elements_new_time_step
        """
        print "Test " + __name__
        mask = self.mask
        mock_compute_dt.return_value = np.zeros(self.nbr_cells)[mask]

        self.my_cells.compute_new_time_step(mask)
        mock_compute_dt.assert_any_call(DataContainer().numeric.cfl,
                                        DataContainer().numeric.cfl_pseudo,
                                        self.my_cells.density.current_value[mask],
                                        self.my_cells.density.new_value[mask],
                                        self.my_cells.size_t_plus_dt[mask],
                                        self.my_cells.sound_velocity.current_value[mask],
                                        self.my_cells.pseudo.current_value[mask],
                                        self.my_cells.pseudo.new_value[mask])
        print "__[OK]"

    def test_compute_complete_stress_tensor(self):
        """
        Test de la méthode compute_complete_stress_tensor
        """
        print "Test " + __name__
        mask = self.mask
        self.my_cells.energy.current_value = np.copy(self.test_variables.energy_old)
        self.my_cells._deviatoric_strain_rate = np.copy(self.test_variables.deviatoric_strain_rate)
        self.my_cells._deviatoric_stress_new = np.copy(self.test_variables.deviatoric_stress_new)
        self.my_cells._deviatoric_stress_current = np.copy(self.test_variables.deviatoric_stress_current)
        self.my_cells.pressure.new_value = np.copy(self.test_variables.pressure_new)
        self.my_cells.pseudo.new_value = np.copy(self.test_variables.pseudo_new)
        self.my_cells.compute_complete_stress_tensor(mask)
        np.testing.assert_allclose(self.my_cells.stress[mask], self.test_variables.stress_new[mask])
        print "__[OK]"

    def test_impose_pressure(self):
        """
        Test de impose_pressure
        """
        print "Test " + __name__
        self.my_cells.pressure.new_value = np.copy(self.test_variables.pressure_new)
        self.my_cells.impose_pressure(0, -1.)
        np.testing.assert_array_equal(self.my_cells.pressure.new_value[1:], self.test_variables.pressure_new[1:])
        np.testing.assert_array_equal(self.my_cells.pressure.new_value[0], np.array([-1.]))
        print "__[OK]"

if __name__ == "__main__":
    unittest.main()
