#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Cell module unit tests
"""
import numpy as np
import unittest
from xvof.src.cell.cell import Cell
from xvof.src.mesh.topology1d import Topology1D
from xvof.src.node import OneDimensionNode
from xvof.src.utilities.testing import captured_output
from xvof.src.data.data_container import DataContainer
from xvof.src.cell.test_cell.test_variables import TestVariables


class CellTest(unittest.TestCase):

    def setUp(self):
        """
        Tests setup
        """
        data_file_path = "//home/marie/PycharmProjects/XVOF/xvof.src/0_UNITTEST/XDATA_hydro.xml"
        self.test_datacontainer = DataContainer(data_file_path)

        self.test_variables = TestVariables(4, 5)
        self.test_variables.define_hydro_variables()

        self.my_cells = Cell(self.test_variables.nb_cells)
        self.my_nodes = OneDimensionNode(5, self.test_variables.node_coord_init, self.test_variables.node_velocity_init)
        self.my_topo = Topology1D(5, 4)

        # Initilisation des champs new_value
        self.my_cells.pressure.new_value = np.copy(self.test_variables.pressure_old)
        self.my_cells.density.new_value = np.copy(self.test_variables.density_old)
        self.my_cells.sound_velocity.new_value = np.copy(self.test_variables.sound_speed_old)
        self.my_cells.energy.new_value = np.copy(self.test_variables.energy_old)

    def tearDown(self):
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        pass

    def test_get_coordinates(self):
        """
        Test of cell.get_coordinates method
        """
        print "Test get_coordinates (cas hydro) :"
        nbr_cell = self.my_cells.number_of_cells
        cell_coord = Cell.get_coordinates(nbr_cell, self.my_topo, self.my_nodes.xt)
        np.testing.assert_allclose(cell_coord, self.test_variables.cell_coord_old)
        print "__[OK]"

    def test_increment_variables(self):
        """
        Test of cell.increment_variables method
        """
        print "Test increment_variables (cas hydro) : "
        self.my_cells.increment_variables()

        # Vérification des résultats
        np.testing.assert_allclose(self.my_cells.pressure.current_value, self.test_variables.pressure_old)
        np.testing.assert_allclose(self.my_cells.density.current_value, self.test_variables.density_old)
        np.testing.assert_allclose(self.my_cells.sound_velocity.current_value, self.test_variables.sound_speed_old)
        np.testing.assert_allclose(self.my_cells.energy.current_value, self.test_variables.energy_old)
        np.testing.assert_allclose(self.my_cells.size_t, self.test_variables.cell_size_old)
        print "__[OK]"

    def test_str(self):
        """
        Teste la réponse de __str__ de la classe Cell
        """
        print "Test __str__"
        message = self.my_cells.__str__()
        answer = "Number of cells: 4"
        self.assertEqual(message.split(), answer.split())
        print "__[OK]"

    def test_print_infos(self):
        """
        Test of cell.print_infos method
        """
        with captured_output() as (out, err):
            self.my_cells.print_infos()
            output = out.getvalue().strip()
            answer = """Cell 
            ==> number of cells = 4
            ==> size at t = [ 0.  0.  0.  0.]
            ==> size at t+dt = [ 0.015  0.01   0.05  0.25]
            ==> density at t = [ 8930.  8930.  8930.  8930.]
            ==> density at t+dt = [ 8500.  4000.  8000.  4000.]
            ==> pressure at t = [ 100000.  100000.  100000.  100000.  100000.]
            ==> pressure at t+dt = [  2.00000000e+09  -5.00000000e+08   1.25000000e+09  0.50000000e+09]
            ==> internal energy at t = [ 6.719465  6.719465  6.719465  6.719465 ]
            ==> internal energy at t+dt = [-1.  0.5   6.  4.]
            ==> sound velocity at t = [ 0.  0.  0.  0.]
            ==> sound velocity at t+dt = [ 200.  300. 250.  300.]"""
            self.assertEqual(output.split(), answer.split())

if __name__ == "__main__":
    unittest.main()
