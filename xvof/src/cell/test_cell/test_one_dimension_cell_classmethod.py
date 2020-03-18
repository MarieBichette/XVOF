#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
one_dimension_cell module unit tests
"""
import numpy as np
import unittest
import os
from xvof.src.cell.one_dimension_cell import OneDimensionCell
from xvof.src.data.data_container import DataContainer
from xvof.src.cell.test_cell.test_variables import TestVariables


class OneDimensionCellTest(unittest.TestCase):

    def setUp(self):
        data_file_path = os.path.realpath(os.path.join(os.getcwd(), "../tests/0_UNITTEST/XDATA_elasto.xml"))
        self.test_datacontainer = DataContainer(data_file_path)
        self.nbr_cells = 4
        self.my_cells = OneDimensionCell(self.nbr_cells)
        self.my_cells.cell_in_target = np.ones(self.nbr_cells, dtype='bool')

        self.test_variables = TestVariables(4, 5)
        self.test_variables.define_elasto_variables()

        # Initilisation des champs new_value
        self.my_cells.pressure.new_value = self.test_variables.pressure_old
        self.my_cells.density.new_value = self.test_variables.density_old
        self.my_cells.sound_velocity.new_value = self.test_variables.sound_speed_old
        self.my_cells.energy.new_value = self.test_variables.energy_old

    def tearDown(self):
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        pass

    def test_apply_equation_of_state(self):
        """
        Test of apply_equation_of_state class method
        """
        print "Test apply_equation_of_state"
        density_current = np.copy(self.test_variables.density_old)
        pressure_current = np.copy(self.test_variables.pressure_old)
        energy_current = np.copy(self.test_variables.energy_old)

        density_new = np.copy(self.test_variables.density_new)
        cson_new = np.zeros([self.nbr_cells])
        pressure_new = np.zeros([self.nbr_cells])
        energy_new = np.zeros([self.nbr_cells])
        pseudo = np.zeros([self.nbr_cells])

        eos = DataContainer().material_target.constitutive_model.eos

        energy_new_value, pressure_new_value, sound_velocity_new_value = \
            OneDimensionCell.apply_equation_of_state(self.my_cells, eos, density_current, density_new,
                                         pressure_current, pressure_new,
                                         energy_current, energy_new, pseudo, cson_new)

        np.testing.assert_allclose(energy_new_value, self.test_variables.energy_new)
        np.testing.assert_allclose(pressure_new_value, self.test_variables.pressure_new)
        np.testing.assert_allclose(sound_velocity_new_value, self.test_variables.sound_speed_new)
        print "__[OK]"

    def test_add_elastic_energy_method(self):
        """
        Test de la méthode add_elastic_energy_method
        """
        print "Test add_elastic_energy_method"
        dt = self.test_variables.dt
        density_current = np.copy(self.test_variables.density_old)
        density_new = np.copy(self.test_variables.density_new)
        stress_dev_current = np.copy(self.test_variables.deviatoric_stress_old)
        stress_dev_new = np.copy(self.test_variables.deviatoric_stress_new)
        strain_rate_dev = np.copy(self.test_variables.strain_rate_dev_new)
        energy_new = OneDimensionCell.add_elastic_energy_method(dt, density_current, density_new,
                                                                stress_dev_current, stress_dev_new, strain_rate_dev)
        np.testing.assert_allclose(energy_new, self.test_variables.energy_new)
        print "__[OK]"

    def test_general_method_deviator_strain_rate(self):
        """
        Test of general_method_deviator_strain_rate
        """
        print "Test general_method_deviator_strain_rate"
        mask = np.ones([self.nbr_cells], dtype=np.bool)
        mask[0] = False
        mask[1] = False
        dt = self.test_variables.dt
        x_new = np.copy(self.test_variables.node_coord_new)
        u_new = np.copy(self.test_variables.node_velocity_new)
        # Reconstruction des array donnés par la topologie
        position_new = np.array([[x_new[0], x_new[1]],
                                 [x_new[1], x_new[2]],
                                 [x_new[2], x_new[3]],
                                 [x_new[3], x_new[4]]]).reshape((4,2))
        vitesse_new = np.array([[u_new[0], u_new[1]],
                                [u_new[1], u_new[2]],
                                [u_new[2], u_new[3]],
                                [u_new[3], u_new[4]]]).reshape((4,2))
        dev_strain_rate = np.zeros(self.test_variables.nb_cells)
        dev_strain_rate[mask] = OneDimensionCell.general_method_deviator_strain_rate(mask, dt, position_new, vitesse_new)
        np.testing.assert_allclose(dev_strain_rate[mask], self.test_variables.strain_rate_dev_new[mask], rtol=1.e-05)
        np.testing.assert_allclose(dev_strain_rate[~mask], np.zeros(2))
        print "__[OK]"

    def test_compute_pseudo(self):
        """
        Test of compute_pseudo class method
        """
        print "Test " + __name__
        rho_new = np.copy(self.test_variables.density_new)
        rho_old = np.copy(self.test_variables.density_old)
        new_size = np.copy(self.test_variables.cell_size_new)
        sound_speed = np.copy(self.test_variables.pseudo_new)
        dt = self.test_variables.dt
        pseudo_a, pseudo_b = 1.5, 0.2
        result = OneDimensionCell.compute_pseudo(dt, rho_old, rho_new, new_size, sound_speed, pseudo_a, pseudo_b)
        np.testing.assert_allclose(result, self.test_variables.pseudo_new)
        print "__[OK]"

    def test_compute_time_step(self):
        """
        Test of compute_time_step class method
        """
        print "Test " + __name__
        cfl, cfl_pseudo = 0.25, 0.1
        rho_new = np.copy(self.test_variables.density_new)
        rho_old = np.copy(self.test_variables.density_old)
        new_size = np.copy(self.test_variables.cell_size_new)
        sound_speed = np.copy(self.test_variables.pseudo_new)
        pseudo_old = np.copy(self.test_variables.pseudo_old)
        pseudo_new = np.copy(self.test_variables.pseudo_new)
        result = OneDimensionCell.compute_time_step(cfl, cfl_pseudo, rho_old, rho_new, new_size, sound_speed,
                                                    pseudo_old, pseudo_new)
        np.testing.assert_allclose(result, self.test_variables.time_step_new)
        print "__[OK]"

if __name__ == "__main__":
    unittest.main()
