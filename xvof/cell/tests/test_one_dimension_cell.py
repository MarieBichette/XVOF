#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
one_dimension_cell module unit tests
"""
import numpy as np
import unittest
import mock

from xvof.cell.one_dimension_cell import OneDimensionCell as Cell


class OneDimensionCellTest(unittest.TestCase):

    def setUp(self):
        """Test setup"""
        self.topology = mock.Mock(name="TopologyMock")
        self.topology.dimension = 2
        self.topology.nodes_belonging_to_cell = np.array([[0, 1], [1, 2], [2, 3]], dtype=np.int64)
        self.my_cell = Cell(3)

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

    def test_compute_size(self):
        """
        Test of compute_size method
        """
        self.my_cell.compute_size(self.topology, np.array([-0.5, 0.1, 0.2, 0.35])) 
        np.testing.assert_allclose(self.my_cell.size_t, [0.6, 0.1, 0.15])
    
    def test_compute_new_size(self):
        """
        Test of compute_new_size method
        """
        self.my_cell.compute_new_size(self.topology, np.array([-0.5, 0.1, 0.2, 0.35]), np.array([True, False, True])) 
        np.testing.assert_allclose(self.my_cell.size_t_plus_dt, [0.6, 0., 0.15])

    def test_compute_mass(self):
        """
        Test of compute_mass method
        """
        self.my_cell.compute_size(self.topology, np.array([-0.5, 0.1, 0.2, 0.35])) 
        self.my_cell.compute_mass()
        np.testing.assert_allclose(self.my_cell.mass, [1.5322804 , 0.25538007, 0.3830701])

    def test_compute_new_density(self):
        """
        Test of compute_new_density method 
        """
        self.my_cell.compute_size(self.topology, np.array([-0.5, 0.1, 0.2, 0.35])) 
        self.my_cell.compute_new_size(self.topology, np.array([-0.6, 0.2, 0.1, 0.37]), np.array([True, True, True]))
        self.my_cell.compute_new_density(np.array([True, True, True]))
        np.testing.assert_allclose(self.my_cell.density.new_value, np.array([6096.75, 8129., 4516.11111111]))



if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
