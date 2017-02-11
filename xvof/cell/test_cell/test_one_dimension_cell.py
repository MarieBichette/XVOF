#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
one_dimension_cell module unit tests
"""
import numpy as np
import unittest
import mock

from xvof.cell.one_dimension_cell import OneDimensionCell as Cell
from xvof.mesh.topology1d import Topology1D
from xvof.data.data_container import geometrical_props
from xvof.fields.field import Field
from xvof.data.data_container import DataContainer


class OneDimensionCellTest(unittest.TestCase):

    def setUp(self):
        self.test_cell = Cell(3)
        self.test_datacontainer = DataContainer(self.test_cell._data_path_file)

    def tearDown(self):
        pass

    def test_compute_pseudo(self):
        """
        Test of compute_pseudo class method
        """
        rho_new = np.array([8500., 3500, 2175])
        rho_old = np.array([8700., 3200, 2171])
        new_size = np.array([0.025, 0.01, 0.005])
        sound_speed = np.array([4400, 3200, 1140])
        dt = 1.2e-08
        pseudo_a, pseudo_b = 1.2, 0.25
        result = Cell.compute_pseudo(dt, rho_old, rho_new, new_size, sound_speed, pseudo_a, pseudo_b) 
        np.testing.assert_allclose(result, [0.00000000e+00, 2.25427729e+13, 2.00897590e+09])

    def test_compute_time_step(self):
        """
        Test of compute_time_step class method
        """
        cfl, cfl_pseudo = 0.25, 0.1
        rho_new = np.array([8500., 2175., 3500.])
        rho_old = np.array([8700., 2174.9, 3200])
        new_size = np.array([0.025, 0.01, 0.005])
        sound_speed = np.array([4400., 3200., 1140.])
        pseudo_old = np.array([1.0e+09, 0.5e+08, 0.3e+08])
        pseudo_new = np.array([1.5e+09, 1.5e+08, 0.])
        result = Cell.compute_time_step(cfl, cfl_pseudo, rho_old, rho_new, new_size, sound_speed, pseudo_old, pseudo_new)    
        np.testing.assert_allclose(result, [1.41137110e-06, 7.81250000e-07, 1.09649123e-06])

    @mock.patch.object(Topology1D, "nodes_belonging_to_cell", new_callable=mock.PropertyMock, return_value=np.array([[0, 1], [1, 2], [2, 3]]))  
    def test_compute_size(self, mock_connectivity):
        """
        Test of compute_size method
        """
        self.test_cell.compute_size(Topology1D(4, 3), np.array([-0.5, 0.1, 0.2, 0.35])) 
        np.testing.assert_allclose(self.test_cell.size_t, [0.6, 0.1, 0.15])
   
    @mock.patch.object(Topology1D, "nodes_belonging_to_cell", new_callable=mock.PropertyMock, return_value=np.array([[0, 1], [1, 2], [2, 3]]))  
    def test_compute_new_size(self, mock_connectivity):
        """
        Test of compute_new_size method
        """
        self.test_cell.compute_new_size(Topology1D(4, 3), np.array([-0.5, 0.1, 0.2, 0.35]), np.array([True, False, True])) 
        np.testing.assert_allclose(self.test_cell.size_t_plus_dt, [0.6, 0., 0.15])

    @mock.patch.object(Cell, "size_t", new_callable=mock.PropertyMock, return_value=np.array([0.6, 0.1, 0.15])) 
    @mock.patch.object(geometrical_props, "section", new_callable=mock.PropertyMock, return_value=0.0003141592653589793) 
    @mock.patch.object(Field, "current_value", new_callable=mock.PropertyMock, return_value=np.array([8129., 8129., 8129.]))
    def test_compute_mass(self, mock_density_filed, mock_geom_props, mock_cell):
        """
        Test of compute_mass method
        """
        self.test_cell.compute_mass()
        np.testing.assert_allclose(self.test_cell.mass, [1.5322804 , 0.25538007, 0.3830701])

    @mock.patch.object(Cell, "size_t", new_callable=mock.PropertyMock, return_value=np.array([0.6, 0.1, 0.15]))
    @mock.patch.object(Cell, "size_t_plus_dt", new_callable=mock.PropertyMock, return_value=np.array([0.8, 0.1, 0.27]))
    def test_compute_new_density(self, mock_new_size, mock_size):
        """
        Test of compute_new_density method 
        """
        self.test_cell.compute_new_density(np.array([True, True, True]))
        np.testing.assert_allclose(self.test_cell.density.new_value, np.array([6096.75, 8129., 4516.11111111]))

    @mock.patch.object(Cell, "energy", new_callable=mock.PropertyMock)
    @mock.patch.object(Cell, "pressure", new_callable=mock.PropertyMock)
    @mock.patch.object(Cell, "density", new_callable=mock.PropertyMock)
    @mock.patch.object(Cell, "pseudo", new_callable=mock.PropertyMock)
    @mock.patch.object(Cell, "sound_velocity", new_callable=mock.PropertyMock)
    def test_compute_new_pressure_internal(self, sound_velocity_mock, pseudo_mock, density_mock, pressure_mock, energy_mock):
        """
        Test of compute_new_pressure method with internal solver
        """
        self.test_cell._external_library = None
        self.test_cell.energy.current_value = np.array([1.e+06, 0.5e+05, 2.4e+07])
        self.test_cell.pressure.current_value = np.array([1.5e+09, 0.5e+08, 1.2e+10])
        self.test_cell.density.current_value = np.array([8000., 8500., 9500.])
        self.test_cell.pseudo.current_value = np.array([1.e+08, 0., 2.4e+09])
        self.test_cell.density.new_value = np.array([8120., 8440., 9620.])
        self.test_cell.energy.new_value = np.array([0., 0., 0.])
        self.test_cell.pressure.new_value = np.array([0., 0., 0.])
        self.test_cell.sound_velocity.new_value = np.array([0., 0., 0.])
        self.test_cell.compute_new_pressure(np.array([True, True, True]))
        # Function to vanish
        delta_v = 1. / self.test_cell.density.new_value - 1. / self.test_cell.density.current_value
        func = (self.test_cell.energy.new_value + self.test_cell.pressure.new_value * delta_v / 2. +
                (self.test_cell.pressure.current_value + 2. * self.test_cell.pseudo.current_value) * delta_v / 2.
                - self.test_cell.energy.current_value)
        np.testing.assert_allclose(func, [0., 0., 0.])

    @mock.patch.object(Cell, "energy", new_callable=mock.PropertyMock)
    @mock.patch.object(Cell, "pressure", new_callable=mock.PropertyMock)
    @mock.patch.object(Cell, "density", new_callable=mock.PropertyMock)
    @mock.patch.object(Cell, "pseudo", new_callable=mock.PropertyMock)
    @mock.patch.object(Cell, "sound_velocity", new_callable=mock.PropertyMock)
    def test_compute_new_pressure_external(self, sound_velocity_mock, pseudo_mock, density_mock, pressure_mock, energy_mock):
        """
        Test of compute_new_pressure method with external solver
        """
        self.test_cell._external_library = ""
        self.test_cell.energy.current_value = np.array([1.e+06, 0.5e+05, 2.4e+07])
        self.test_cell.pressure.current_value = np.array([1.5e+09, 0.5e+08, 1.2e+10])
        self.test_cell.density.current_value = np.array([8000., 8500., 9500.])
        self.test_cell.pseudo.current_value = np.array([1.e+08, 0., 2.4e+09])
        self.test_cell.density.new_value = np.array([8120., 8440., 9620.])
        self.test_cell.energy.new_value = np.array([0., 0., 0.])
        self.test_cell.pressure.new_value = np.array([0., 0., 0.])
        self.test_cell.sound_velocity.new_value = np.array([0., 0., 0.])
        self.test_cell.compute_new_pressure(np.array([True, True, True]))
        # Function to vanish
        delta_v = 1. / self.test_cell.density.new_value - 1. / self.test_cell.density.current_value
        func = (self.test_cell.energy.new_value + self.test_cell.pressure.new_value * delta_v / 2. +
                (self.test_cell.pressure.current_value + 2. * self.test_cell.pseudo.current_value) * delta_v / 2.
                - self.test_cell.energy.current_value)
        np.testing.assert_allclose(func, [0., 0., 0.])

    @mock.patch.object(Cell, "energy", new_callable=mock.PropertyMock)
    @mock.patch.object(Cell, "pressure", new_callable=mock.PropertyMock)
    @mock.patch.object(Cell, "density", new_callable=mock.PropertyMock)
    @mock.patch.object(Cell, "pseudo", new_callable=mock.PropertyMock)
    @mock.patch.object(Cell, "sound_velocity", new_callable=mock.PropertyMock)
    @mock.patch.object(DataContainer, "numeric", new_callable= mock.PropertyMock)
    def test_compute_new_pseudo(self):
        """   Test of compute_new_pseudo     """
        mask = [True, True, True]
        self.test_cell.density.current_value = np.array([8000., 8500., 9500.])
        self.test_cell.density.new_value = np.array([8120., 8440., 9620.])
        self.test_cell.sound_velocity.current_value = np.array([0., 0., 0.])
        self.test_datacontainer.numeric.a_pseudo = 0.2
        self.test_datacontainer.numeric.b_pseudo = 0.0
        self.test_cell.pseudo.new_value[mask] = Cell.compute_pseudo(1.0e+06, self.test_cell.density.current_value[mask],
                                                                      self.test_cell.density.new_value[mask],
                                                                      self.test_cell.size_t_plus_dt[mask],
                                                                      self.test_cell.sound_velocity.current_value[mask],
                                                                      self.test_datacontainer.numeric.a_pseudo,
                                                                      self.test_datacontainer.numeric.b_pseudo)
        np.testing.assert_allclose(self.test_cell.pseudo.new_value[mask], np.array([0., 0., 0.]))


    @mock.patch.object(Cell, "density", new_callable=mock.PropertyMock)
    @mock.patch.object(Cell, "pseudo", new_callable=mock.PropertyMock)
    @mock.patch.object(Cell, "sound_velocity", new_callable=mock.PropertyMock)
    @mock.patch.object(DataContainer, "numeric", new_callable= mock.PropertyMock)
    @mock.patch.object(Cell, "size_t_plus_dt", new_callable = mock.PropertyMock, return_value = np.array([1.]))
    def test_compute_new_time_step(self, mock_sizetpdt):
        """Test de compute_new_time_step"""
        self.test_cell.density.current_value = np.array([8000., 8500., 9500.])
        self.test_cell.density.new_value = np.array([8120., 8440., 9620.])
        self.test_cell.sound_velocity.new_value = np.array([0., 0., 0.])
        self.test_cell.pseudo.current_value = np.array([0., 0., 0.])
        self.test_cell.pseudo.new_value = np.array([0., 0., 0.])
        cfl = self.test_datacontainer.numeric.cfl
        cfl_pseudo = self.test_datacontainer.numeric.cfl_pseudo
        mask = np.array([False, True, False])
        self.test_cell.compute_new_time_step(mask)
        # compute_time_step(cfl, cfl_pseudo, self.test_cell.density.current_value, self.test_cell.density.new_value,
        #                                    self.test_cell.size_t_plus_dt, self.test_cell.sound_velocity.new_value,
        #                                    self.test_cell.pseudo.current_value, self.test_cell.pseudo.new_value)
        np.testing.assert_array_equal(self.test_cell._dt, 4.) #calcul à faire

    @mock.patch.object(Cell, "pressure", new_callable = mock.PropertyMock)
    def test_impose_pressure(self, mock_pressure):
        """Test de impose_pressure"""
        self.test_cell.pressure.new_value = np.array([1., 2., 3.])
        self.test_cell.impose_pressure(1, 4.)
        np.testing.assert_array_equal(self.test_cell.pressure.new_value, np.array([1., 4., 3.]))


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
