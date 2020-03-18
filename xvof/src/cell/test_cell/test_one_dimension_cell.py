#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
one_dimension_cell module unit tests
"""
import numpy as np
import unittest
import mock
import os

from xvof.src.cell.one_dimension_cell import OneDimensionCell as Cell
from xvof.src.mesh.topology1d import Topology1D
from xvof.src.data.data_container import geometrical_props
from xvof.src.fields.field import Field
from xvof.src.data.data_container import DataContainer


class OneDimensionCellTest(unittest.TestCase):

    def setUp(self):
        self.test_cell = Cell(3)
        data_file_path = os.path.realpath(os.path.join(os.path.dirname(__file__), "../../../tests/0_UNITTEST/XDATA.xml"))
        self.test_datacontainer = DataContainer(data_file_path)

    def tearDown(self):
        pass

    @mock.patch.object(Cell, "size_t", new_callable=mock.PropertyMock, return_value=np.array([0.6, 0.1, 0.15])) 
    @mock.patch.object(geometrical_props, "section", new_callable=mock.PropertyMock, return_value=0.0003141592653589793) 
    @mock.patch.object(Field, "current_value", new_callable=mock.PropertyMock, return_value=np.array([8129., 8129., 8129.]))
    def test_compute_mass(self, mock_density_filed, mock_geom_props, mock_cell):
        """
        Test of compute_mass method
        """
        self.test_cell.compute_mass()
        np.testing.assert_allclose(self.test_cell.mass, [1.5322804 , 0.25538007, 0.3830701])

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

    @mock.patch.object(Cell, "pressure", new_callable = mock.PropertyMock)
    def test_impose_pressure(self, mock_pressure):
        """Test de impose_pressure"""
        self.test_cell.pressure.new_value = np.array([1., 2., 3.])
        self.test_cell.impose_pressure(1, 4.)
        np.testing.assert_array_equal(self.test_cell.pressure.new_value, np.array([1., 4., 3.]))


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
