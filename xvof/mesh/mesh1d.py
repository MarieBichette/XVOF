#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Base class for one dimensional mesh
"""
import numpy as np

from xvof.cell.one_dimension_cell import OneDimensionCell
from xvof.data.data_container import DataContainer
from xvof.mesh.topology1d import Topology1D
from xvof.node.one_dimension_node import OneDimensionNode


class Mesh1d(object):
    """
    This class defines a one dimensional mesh
    """

    def __init__(self, initial_coordinates, initial_velocities):
        if np.shape(initial_coordinates) != np.shape(initial_velocities):
            message = "Initial velocity and coordinates vector doesn't have the same shape!"
            raise ValueError(message)
        if np.shape(initial_coordinates)[1] != 1:
            message = ("""A 1D mesh must have one dimensional vector which is not the case"""
                       """ for initial coordinates vector!""")
            raise ValueError(message)
        # ---------------------------------------------
        # Nodes creation
        # ---------------------------------------------
        nbr_nodes = np.shape(initial_coordinates)[0]
        self.nodes = OneDimensionNode(nbr_nodes, initial_coordinates, initial_velocities,
                                      section=DataContainer().geometric.section)
        # ---------------------------------------------
        # Cells creation
        # ---------------------------------------------
        nbr_cells = nbr_nodes - 1
        self.cells = OneDimensionCell(nbr_cells)
        self._all_cells = np.zeros(self.cells.number_of_cells, dtype=np.bool, order='C')
        self._all_cells[:] = True
        # ---------------------------------------------
        # Topology creation
        # ---------------------------------------------
        self.__topology = Topology1D(nbr_nodes, nbr_cells)
        # ---------------------------------------------
        # Ruptured cells vector
        # ---------------------------------------------
        self.__ruptured_cells = np.zeros(self.cells.number_of_cells, dtype=np.bool, order='C')

    def compute_nodes_masses(self):
        """ Nodal mass computation """
        nb_nodes_per_cell = np.zeros([self.cells.number_of_cells, ], dtype=np.int, order='C')
        nb_nodes_per_cell[:] = 2
        self.nodes.calculer_masse_wilkins(self.__topology, self.cells.mass, nb_nodes_per_cell)

    def compute_new_nodes_velocities(self, delta_t):
        """
        Computation of nodes velocities at t+dt

        :var delta_t: time step
        :type delta_t: float
        """
        self.nodes.compute_new_velocity(delta_t)

    def compute_new_nodes_coordinates(self, delta_t):
        """
        Computation of nodes coordinates at t+dt

        :var delta_t: time step
        :type delta_t: float
        """
        self.nodes.compute_new_coodinates(delta_t)

    def compute_cells_sizes(self):
        """
        Computation of cells sizes at t
        """
        self.cells.compute_size(self.__topology, self.nodes.xt)

    def compute_new_cells_sizes(self, delta_t):
        """
        Computation of cells sizes at t+dt
        """
        self.cells.compute_new_size(self.__topology, self.nodes.xtpdt, self._all_cells, time_step=delta_t)

    def compute_new_cells_densities(self):
        """
        Computation of cells densities at t+dt
        """
        self.cells.compute_new_density(mask=self._all_cells)

    def compute_new_cells_pressures(self):
        """
        Computation of cells pressure at t+dt
        """
        self.cells.compute_new_pressure(mask=~self.__ruptured_cells)

    def compute_new_cells_pseudo_viscosity(self, delta_t):
        """
        Computation of cells artificial viscosity at t+dt

        :var delta_t: time step
        :type delta_t: float
        """
        self.cells.compute_new_pseudo(delta_t, mask=self._all_cells)

    def compute_new_nodes_forces(self):
        """
        Computation of nodes forces at t+dt
        """
        self.nodes.compute_new_force(self.__topology, self.cells.pressure.new_value, self.cells.pseudo.new_value)

    def increment(self):
        """
        Moving to next time step
        """
        self.nodes.increment()
        self.cells.increment_variables()

    def compute_new_time_step(self):
        """
        Computation of new time step
        """
        self.cells.compute_new_time_step(mask=self._all_cells)
        return self.cells.dt.min()

    def apply_pressure(self, surface, pressure):
        """
        Apply a given pressure on left or right boundary

        :var surface: name of the surface where pressure has to be imposed
        :var pressure: value of the pressure to impose
        :type surface: str ('left' | 'right')
        :type pressure: float
        """
        if surface.lower() not in ("left", "right"):
            raise (ValueError("One dimensional mesh : only 'left' or 'right' boundaries are possibles!"))
        if surface.lower() == 'left':
            self.nodes.apply_pressure(0, pressure)
        else:
            self.nodes.apply_pressure(-1, -pressure)

    def get_ruptured_cells(self, rupture_criterion):
        """
        Find the cells where the rupture criterion is checked and store them

        :var rupture_criterion: rupture criterion
        :type rupture_criterion: RuptureCriterion
        """
        self.__ruptured_cells = np.logical_or(self.__ruptured_cells, rupture_criterion.checkCriterion(self.cells))

    def apply_rupture_treatment(self, treatment):
        """
        Apply the rupture treatment on the cells enforcing the rupture criterion

        :var treatment: rupture treatment
        :type treatment: RuptureTreatment
        """
        treatment.applyTreatment(self.cells, self.__ruptured_cells)

    @property
    def velocity_field(self):
        """
        Node velocity field
        """
        return self.nodes.upundemi

    @property
    def nodes_coordinates(self):
        """
        Nodes coordinates
        """
        return self.nodes.xtpdt

    @property
    def cells_coordinates(self):
        """
        Cells coordinates (coordinates of cells centers)
        """
        return self.cells.get_coordinates(self.cells.number_of_cells, self.__topology, self.nodes.xt)

    @property
    def pressure_field(self):
        """
        Pressure field
        """
        return self.cells.pressure_field

    @property
    def density_field(self):
        """
        Density field
        """
        return self.cells.density_field

    @property
    def energy_field(self):
        """
        Internal energy field
        """
        return self.cells.energy_field

    @property
    def artificial_viscosity_field(self):
        """
        Artificial viscosity field
        """
        return self.cells.artificial_viscosity_field
