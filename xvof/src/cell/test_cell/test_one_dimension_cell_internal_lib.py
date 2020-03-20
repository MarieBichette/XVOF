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


# TODO : move to the package solver !

class OneDimensionCellInternalLibTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("Appel de setUpForClass")
        data_file_path = os.path.join(os.path.dirname(__file__), "../../../tests/0_UNITTEST/XDATA_hydro.xml")
        DataContainer(data_file_path)
        print(DataContainer().material_target.constitutive_model.elasticity_model is not None)
        print(DataContainer().material_projectile.constitutive_model.elasticity_model is not None)
        print(DataContainer().material_target.constitutive_model.plasticity_model is not None)
        print(DataContainer().material_projectile.constitutive_model.plasticity_model is not None)

    @classmethod
    def tearDownClass(cls):
        print("Appel de tearDownForClass")
        DataContainer.clear()
        print("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")

        pass

    def setUp(self):
        self.test_cell = Cell(3)

    def tearDown(self):
        pass

    def test_compute_new_pressure_internal(self):
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
        self.test_cell.compute_new_pressure(np.array([True, True, True]), 1.e-6)
        # Function to vanish
        delta_v = 1. / self.test_cell.density.new_value - 1. / self.test_cell.density.current_value
        func = (self.test_cell.energy.new_value + self.test_cell.pressure.new_value * delta_v / 2. +
                (self.test_cell.pressure.current_value + 2. * self.test_cell.pseudo.current_value) * delta_v / 2.
                - self.test_cell.energy.current_value)
        np.testing.assert_allclose(func, np.array([-1001570.197044, -49979.091163, -24011029.653135]))


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
