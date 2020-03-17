#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Cell module unit tests
"""
import numpy as np
import unittest
import mock

from xvof.src.cell.cell import Cell
from xvof.src.mesh.topology1d import Topology1D
from xvof.src.node import OneDimensionNode
from xvof.src.utilities.testing import captured_output


class CellTest(unittest.TestCase):

    def setUp(self):
        """
        Tests setup
        """
        self.my_cell = Cell(3)

    def tearDown(self):
        pass

    def test_get_coordinates(self):
        """
        Test of cell.get_coordinates method
        """
        nodes = OneDimensionNode(4, np.array([[-0.5], [0.1], [0.2], [0.35]]), np.array([[0.], [0.], [0.], [0.]]))
        my_topo = Topology1D(4, 3)
        cell_coord = Cell.get_coordinates(3, my_topo, nodes.xt)
        np.testing.assert_allclose(cell_coord, np.array([[-0.2], [0.15], [0.275]]))


    def test_increment_variables(self):
        """
        Test of cell.increment_variables method
        """
        self.my_cell.pressure.new_value = np.array([2.0e+09, -0.5e+09, 1.25e+09])
        self.my_cell.density.new_value = np.array([8500., 4500, 9500.])
        self.my_cell.sound_velocity.new_value = np.array([440., 210, -110])
        self.my_cell.energy.new_value = np.array([-1.0e+06, 0.5e+05, 1.5e+05]) 
        self.my_cell._size_t_plus_dt = np.array([0.015, 0.01, 0.05])
        self.my_cell.increment_variables()
        np.testing.assert_allclose(self.my_cell.pressure.current_value, np.array([2.0e+09, -0.5e+09, 1.25e+09]) )
        np.testing.assert_allclose(self.my_cell.density.current_value, np.array([8500., 4500, 9500.]))
        np.testing.assert_allclose(self.my_cell.sound_velocity.current_value, np.array([440., 210, -110]))
        np.testing.assert_allclose(self.my_cell.energy.current_value, np.array([-1.0e+06, 0.5e+05, 1.5e+05]) )
        np.testing.assert_allclose(self.my_cell.size_t, np.array([0.015, 0.01, 0.05]))

    def test_str(self):
        """Teste la réponse de __str__ de la classe Cell"""
        message = self.my_cell.__str__()
        answer = "Number of cells: 3"
        self.assertEqual(message.split(), answer.split())

    def test_print_infos(self):
        """
        Test of cell.print_infos method
        """
        self.my_cell.pressure.new_value = np.array([2.0e+09, -0.5e+09, 1.25e+09])
        self.my_cell.density.new_value = np.array([8500., 4500, 9500.])
        self.my_cell.sound_velocity.new_value = np.array([440., 210, -110])
        self.my_cell.energy.new_value = np.array([-1.0e+06, 0.5e+05, 1.5e+05]) 
        self.my_cell._size_t_plus_dt = np.array([0.015, 0.01, 0.05])
        with captured_output() as (out, err):
            self.my_cell.print_infos()
            output = out.getvalue().strip()
            answer = """Cell 
            ==> number of cells = 3
            ==> size at t = [ 0.  0.  0.]
            ==> size at t+dt = [ 0.015  0.01   0.05 ]
            ==> density at t = [ 8129.  8129.  8129.]
            ==> density at t+dt = [ 8500.  4500.  9500.]
            ==> pressure at t = [ 100006.2096  100006.2096  100006.2096]
            ==> pressure at t+dt = [  2.00000000e+09  -5.00000000e+08   1.25000000e+09]
            ==> internal energy at t = [ 7.689  7.689  7.689]
            ==> internal energy at t+dt = [-1000000.    50000.   150000.]
            ==> sound velocity at t = [ 0.  0.  0.]
            ==> sound velocity at t+dt = [ 440.  210. -110.]"""
            self.assertEqual(output.split(), answer.split())

if __name__ == "__main__":
    unittest.main()
