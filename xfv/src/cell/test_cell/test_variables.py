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


class TestVariables:

    def __init__(self, nb_cells, nb_nodes):
        """

        :param nb_cells:
        :param nbr_nodes:
        """
        self.dt = 1.
        self.nb_cells = nb_cells
        self.nb_nodes = nb_nodes

        self.density_init = np.ones(self.nb_cells) * 8930.
        self.pressure_init = np.ones(self.nb_cells) * 1.e+05
        self.energy_init = np.ones(self.nb_cells) * 6.719465
        self.pseudo_init = np.ones(self.nb_cells) * 0.
        self.sound_speed_init = np.ones(self.nb_cells) * 8930.
        self.shear_modulus_init = np.ones(self.nb_cells) * 47.7e+09
        self.yield_stress_init = np.ones(self.nb_cells) * 120.e+06

        self.deviatoric_stress_init = np.array([2000., 4000, 8000., 4000.])
        self.stress_init = np.array([2000., 4000, 8000., 4000.])
        self.stress_xx_init = np.array([2000., 4000, 8000., 4000.])
        self.strain_rate_dev_init = np.array([2000., 4000, 8000., 4000.])

        self.node_coord_init = np.array([[-0.5, ], [0.1, ], [0.2, ], [0.35, ], [0.65, ]])
        self.node_coord_old = np.array([[-0.35, ], [0.15, ], [0.2, ], [0.25, ], [0.5, ]])
        self.node_coord_new = np.array([[-0.25, ], [0.1, ], [0.2, ], [0.45, ], [0.85, ]])

        self.node_velocity_init = np.zeros([self.nb_nodes, 1])
        self.node_velocity_old = np.array([[0.15, ], [0.05, ], [0., ], [-0.1, ], [-0.15, ]])
        self.node_velocity_new = np.array([[0.1, ], [-0.05, ], [0., ], [0.2, ], [0.3, ]])

        self.cell_coord_init = np.zeros([self.nb_nodes, 1])
        self.cell_coord_old = np.array([[0.15, ], [0.05, ], [0., ], [-0.1, ], [-0.15, ]])
        self.cell_coord_new = np.array([[0.1, ], [-0.05, ], [0., ], [0.2, ], [0.3, ]])

        self.time_step_new = np.array([2000., 4000, 8000., 4000.])
        self.cell_mass = np.array([5358., 893, 1339.5, 2679.])

    def define_hydro_variables(self):
        """
        Définition des grandeurs old et new pour un cas hydro
        :return:
        """
        self.density_old = np.array([10716., 17860, 26790., 10716.])
        self.density_new = np.array([15308.5143, 8930, 5358., 6697.5])

        self.pressure_old = np.array([2000., 4000, 8000., 4000.])
        self.pressure_new = np.array([2000., 4000, 8000., 4000.])

        self.energy_old = np.array([2000., 4000, 8000., 4000.])
        self.energy_new = np.array([2000., 4000, 8000., 4000.])

        self.pseudo_old = np.array([2000., 4000, 8000., 4000.])
        self.pseudo_new = np.array([2000., 4000, 8000., 4000.])

        self.sound_speed_old = np.array([2000., 4000, 8000., 4000.])
        self.sound_speed_new = np.array([2000., 4000, 8000., 4000.])

        self.deviatoric_stress_old = np.zeros([self.nb_cells, 3])
        self.deviatoric_stress_new = np.zeros([self.nb_cells, 3])
        self.deviatoric_strain_rate = np.zeros([self.nb_cells, 3])

        self.stress_old = np.array(
            [[-(self.pseudo_old[0] + self.pressure_old[0]), -(self.pseudo_old[0] + self.pressure_old[0]), -(self.pseudo_old[0] + self.pressure_old[0])],
             [-(self.pseudo_old[1] + self.pressure_old[1]), -(self.pseudo_old[1] + self.pressure_old[1]), -(self.pseudo_old[1] + self.pressure_old[1])],
             [-(self.pseudo_old[2] + self.pressure_old[2]), -(self.pseudo_old[2] + self.pressure_old[2]), -(self.pseudo_old[2] + self.pressure_old[2])],
             [-(self.pseudo_old[3] + self.pressure_old[3]), -(self.pseudo_old[3] + self.pressure_old[3]), -(self.pseudo_old[3] + self.pressure_old[3])]])

        self.stress_new = np.array(
            [[-(self.pseudo_new[0] + self.pressure_new[0]), -(self.pseudo_new[0] + self.pressure_new[0]),
              -(self.pseudo_new[0] + self.pressure_new[0])],
             [-(self.pseudo_new[1] + self.pressure_new[1]), -(self.pseudo_new[1] + self.pressure_new[1]),
              -(self.pseudo_new[1] + self.pressure_new[1])],
             [-(self.pseudo_new[2] + self.pressure_new[2]), -(self.pseudo_new[2] + self.pressure_new[2]),
              -(self.pseudo_new[2] + self.pressure_new[2])],
             [-(self.pseudo_new[3] + self.pressure_new[3]), -(self.pseudo_new[3] + self.pressure_new[3]),
              -(self.pseudo_new[3] + self.pressure_new[3])]])

        self.cell_size_init = np.array([0.6, 0.1, 0.15, 0.3])
        self.cell_size_old = np.array([0.5, 0.05, 0.05, 0.25])
        self.cell_size_new = np.array([0.35, 0.1, 0.25, 0.4])

    def define_elasto_variables(self):
        """
        Définition des grandeurs old et new pour un cas elasto
        :return:
        """
        self.density_old = np.array([2000., 4000, 8000., 4000.])
        self.density_new = np.array([2000., 4000, 8000., 4000.])

        self.pressure_old = np.array([2000., 4000, 8000., 4000.])
        self.pressure_new = np.array([2000., 4000, 8000., 4000.])

        self.energy_old = np.array([2000., 4000, 8000., 4000.])
        self.energy_new = np.array([2000., 4000, 8000., 4000.])

        self.pseudo_old = np.array([2000., 4000, 8000., 4000.])
        self.pseudo_new = np.array([2000., 4000, 8000., 4000.])

        self.sound_speed_old = np.array([2000., 4000, 8000., 4000.])
        self.sound_speed_new = np.array([2000., 4000, 8000., 4000.])

        self.shear_modulus_old = np.array([2000., 4000, 8000., 4000.])
        self.shear_modulus_new = np.array([2000., 4000, 8000., 4000.])

        self.yield_stress_old = np.array([2000., 4000, 8000., 4000.])
        self.yield_stress_new = np.array([2000., 4000, 8000., 4000.])

        self.deviatoric_stress_old = np.array([[2000., 4000, 8000.], [2000., 4000, 8000.], [2000., 4000, 8000.], [2000., 4000, 8000.]])
        self.deviatoric_stress_new = np.array([[2000., 4000, 8000.], [2000., 4000, 8000.], [2000., 4000, 8000.], [2000., 4000, 8000.]])

        self.strain_rate_dev_old = np.array([[2000., 4000, 8000.], [2000., 4000, 8000.], [2000., 4000, 8000.], [2000., 4000, 8000.]])
        self.strain_rate_dev_new = np.array([[2000., 4000, 8000.], [2000., 4000, 8000.], [2000., 4000, 8000.], [2000., 4000, 8000.]])

        self.stress_old = np.array([[2000., 4000, 8000.], [2000., 4000, 8000.], [2000., 4000, 8000.], [2000., 4000, 8000.]])
        self.stress_new = np.array([[2000., 4000, 8000.], [2000., 4000, 8000.], [2000., 4000, 8000.], [2000., 4000, 8000.]])

        self.stess_xx_old = np.array([2000., 4000, 8000., 4000.])
        self.stess_xx_new = np.array([2000., 4000, 8000., 4000.])

        self.cell_size_old = np.array([0.5, 0.05, 0.05, 0.25])
        self.cell_size_new = np.array([0.35, 0.1, 0.25, 0.4])

    def define_epp_variables(self):
        """
        Définition des grandeurs old et new pour un cas elasto
        :return:
        """
        self.define_elasto_variables()
        self.deviatoric_stress_new_after_plastic_correction = np.array([[2000., 4000, 8000.], [2000., 4000, 8000.], [2000., 4000, 8000.], [2000., 4000, 8000.]])
        self.equivalent_plastic_strain_rate = np.array([2000., 4000, 8000., 4000.])