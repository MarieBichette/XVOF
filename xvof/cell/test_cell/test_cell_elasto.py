#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Cell module unit tests
"""
import numpy as np
import unittest
from xvof.cell.cell import Cell
from xvof.mesh.topology1d import Topology1D
from xvof.node import OneDimensionNode
from xvof.data.data_container import DataContainer
from xvof.cell.test_cell.test_variables import TestVariables


class CellTest(unittest.TestCase):

    def setUp(self):
        """
        Tests setup
        """
        data_file_path = "//home/marie/PycharmProjects/XVOF/xvof/0_UNITTEST/XDATA_elasto.xml"
        self.test_datacontainer = DataContainer(data_file_path)

        self.test_variables = TestVariables(4, 5)
        self.test_variables.define_elasto_variables()

        self.my_cells = Cell(4)
        self.my_nodes = OneDimensionNode(5, self.test_variables.node_coord_init, self.test_variables.node_velocity_init)
        self.my_topo = Topology1D(5, 4)

        # Initilisation des champs new_value
        self.my_cells.pressure.new_value = np.copy(self.test_variables.pressure_old)
        self.my_cells.density.new_value = np.copy(self.test_variables.density_old)
        self.my_cells.sound_velocity.new_value = np.copy(self.test_variables.sound_speed_old)
        self.my_cells.energy.new_value = np.copy(self.test_variables.energy_old)
        self.my_cells.shear_modulus.new_value = np.copy(self.test_variables.shear_modulus_old)
        self.my_cells.yield_stress.new_value = np.copy(self.test_variables.yield_stress_old)
        self.my_cells.deviatoric_stress_new = np.copy(self.test_variables.deviatoric_stress_old)

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
        np.testing.assert_allclose(cell_coord, self.test_variables.cell_size_old)
        print "__[OK]"

    def test_increment_variables(self):
        """
        Test of cell.increment_variables method
        """
        print "Test increment_variables (cas hydro) : "
        self.my_cells.increment_variables()

        # Vérification des résultats
        np.testing.assert_allclose(self.my_cells.shear_modulus.current_value, self.test_variables.shear_modulus_old)
        np.testing.assert_allclose(self.my_cells.yield_stress.current_value, self.test_variables.yield_stress_old)
        print "__[OK]"

if __name__ == "__main__":
    unittest.main()
