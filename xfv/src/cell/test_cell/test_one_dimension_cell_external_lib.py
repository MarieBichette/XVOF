#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
one_dimension_cell module unit tests
"""
import os
import unittest
import numpy as np

from xfv.src.cell.one_dimension_cell import OneDimensionCell as Cell
from xfv.src.data.data_container import DataContainer

# TODO : fix external library path

class OneDimensionCellTest(unittest.TestCase):
    """Test the OneDimensionCell class in case of external library use"""
    @unittest.skip("Solve path to external lib")
    def setUp(self):
        data_file_path = os.path.realpath(
            os.path.join(os.path.dirname(__file__), "../../../tests/0_UNITTEST/XDATA.xml"))
        self.test_datacontainer = DataContainer(data_file_path)
        self.test_cell = Cell(3)

    @unittest.skip("Solve path to external lib")
    def tearDown(self):
        DataContainer.clear()

    @unittest.skip("Solve path to external lib")
    def test_compute_new_pressure_external(self):
        """
        Test of compute_new_pressure method with external solver
        """
        self.test_cell._external_library = ""  # pylint:disable=protected-access
        self.test_cell.energy.current_value = np.array([1.e+06, 0.5e+05, 2.4e+07])
        self.test_cell.pressure.current_value = np.array([1.5e+09, 0.5e+08, 1.2e+10])
        self.test_cell.density.current_value = np.array([8000., 8500., 9500.])
        self.test_cell.pseudo.current_value = np.array([1.e+08, 0., 2.4e+09])
        self.test_cell.density.new_value = np.array([8120., 8440., 9620.])
        self.test_cell.energy.new_value = np.array([0., 0., 0.])
        self.test_cell.pressure.new_value = np.array([0., 0., 0.])
        self.test_cell.sound_velocity.new_value = np.array([0., 0., 0.])
        self.test_cell.compute_new_pressure(np.array([True, True, True]), 1.e-6)
        # Function to vanish
        delta_v = 1. / self.test_cell.density.new_value - 1. / self.test_cell.density.current_value
        func = (self.test_cell.energy.new_value + self.test_cell.pressure.new_value * delta_v / 2. +
                (self.test_cell.pressure.current_value + 2. * self.test_cell.pseudo.current_value)
                * delta_v / 2. - self.test_cell.energy.current_value)
        np.testing.assert_allclose(func, [0., 0., 0.])


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
