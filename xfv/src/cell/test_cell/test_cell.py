# -*- coding: utf-8 -*-
# pylint: disable=protected-access
"""
Cell module unit tests
"""
import numpy as np
import unittest
import unittest.mock as mock
import os

from xfv.src.cell.cell import Cell
from xfv.src.mesh.topology1d import Topology1D
from xfv.src.node import OneDimensionNode
from xfv.src.utilities.testing import captured_output
from xfv.src.data.data_container import DataContainer


class CellTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Tests setup for class
        """
        data_file_path = os.path.join(os.path.dirname(__file__),
                                      "../../../tests/0_UNITTEST/XDATA_elasto.json")
        DataContainer(data_file_path)

    @classmethod
    def tearDownClass(cls):
        DataContainer.clear()
        print("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")

    def setUp(self):
        """
        Tests setup
        """
        self.nbr_cells = 3
        self.my_cells = Cell(self.nbr_cells)

    def tearDown(self):
        pass

    def test_get_coordinates(self):
        """
        Test of cell.get_coordinates method
        """
        nodes = OneDimensionNode(4, np.array([[-0.5], [0.1], [0.2], [0.35]]),
                                 np.array([[0.], [0.], [0.], [0.]]))
        my_topo = Topology1D(4, 3)
        cell_coord = Cell.get_coordinates(3, my_topo, nodes.xt)
        np.testing.assert_allclose(cell_coord, np.array([[-0.2], [0.15], [0.275]]))

    def test_increment_variables(self):
        """
        Test of cell.increment_variables method
        """
        self.my_cells.pressure.new_value = np.array([2.0e+09, -0.5e+09, 1.25e+09])
        self.my_cells.density.new_value = np.array([8500., 4500, 9500.])
        self.my_cells.sound_velocity.new_value = np.array([440., 210, -110])
        self.my_cells.energy.new_value = np.array([-1.0e+06, 0.5e+05, 1.5e+05])
        self.my_cells._size_t_plus_dt = np.array([0.015, 0.01, 0.05])
        self.my_cells.shear_modulus.new_value = np.array([2., 3., 4.])
        self.my_cells.yield_stress.new_value = np.array([1., 2., 3.])
        self.my_cells.increment_variables()
        np.testing.assert_allclose(self.my_cells.pressure.current_value,
                                   np.array([2.0e+09, -0.5e+09, 1.25e+09]))
        np.testing.assert_allclose(self.my_cells.density.current_value,
                                   np.array([8500., 4500, 9500.]))
        np.testing.assert_allclose(self.my_cells.sound_velocity.current_value,
                                   np.array([440., 210, -110]))
        np.testing.assert_allclose(self.my_cells.energy.current_value,
                                   np.array([-1.0e+06, 0.5e+05, 1.5e+05]))
        np.testing.assert_allclose(self.my_cells.size_t, np.array([0.015, 0.01, 0.05]))
        np.testing.assert_allclose(self.my_cells.shear_modulus.current_value,
                                   np.array([2., 3., 4.]))
        np.testing.assert_allclose(self.my_cells.yield_stress.current_value, np.array([1., 2., 3.]))

    def test_str(self):
        """
        Teste la r�ponse de __str__ de la classe Cell
        """
        message = self.my_cells.__str__()
        answer = "Number of cells: 3"
        self.assertEqual(message.split(), answer.split())

    def test_print_infos(self):
        """
        Test of cell.print_infos method
        """
        self.my_cells.pressure.current_value = np.ones(self.nbr_cells) * 100006.2096
        self.my_cells.pressure.new_value = np.array([2.0e+09, -0.5e+09, 1.25e+09])
        self.my_cells.density.current_value = np.ones(self.nbr_cells) * 8129.
        self.my_cells.density.new_value = np.array([8500., 4500., 9500.])
        self.my_cells.sound_velocity.new_value = np.array([440., 210, -110])
        self.my_cells.energy.current_value = np.ones(self.nbr_cells) * 7.689
        self.my_cells.energy.new_value = np.array([-1.0e+06, 0.5e+05, 1.5e+05])
        self.my_cells._size_t = np.zeros(self.nbr_cells)
        self.my_cells._size_t_plus_dt = np.array([0.015, 0.01, 0.05])
        with captured_output() as (out, err):
            self.my_cells.print_infos()
            output = out.getvalue().strip()
        answer = ("Cell " + os.linesep +
                  "==> number of cells = 3" + os.linesep +
                  "==> size at t = [0. 0. 0.]" + os.linesep +
                  "==> size at t+dt = [0.015 0.01  0.05 ]" + os.linesep +
                  "==> density at t = [8129. 8129. 8129.]" + os.linesep +
                  "==> density at t+dt = [8500. 4500. 9500.]" + os.linesep +
                  "==> pressure at t = [100006.2096 100006.2096 100006.2096]" + os.linesep +
                  "==> pressure at t+dt = [ 2.00e+09 -5.00e+08  1.25e+09]" + os.linesep +
                  "==> internal energy at t = [7.689 7.689 7.689]" + os.linesep +
                  "==> internal energy at t+dt = [-1000000.    50000.   150000.]" + os.linesep +
                  "==> sound velocity at t = [0. 0. 0.]" + os.linesep +
                  "==> sound velocity at t+dt = [ 440.  210. -110.]")
        self.assertEqual(output.split('\n'), answer.split('\n'))


if __name__ == "__main__":
    unittest.main()
