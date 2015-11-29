#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Base class for one dimensional mesh
"""
import numpy as np
from xvof.element.element1d import Element1d
from xvof.mesh.topology1d import Topology1D
from xvof.node.node1d import Node1d
from xvof.utilities.profilingperso import timeit_file


class Mesh1d(object):
    """
    This class defines a one dimensional mesh
    """
    def __init__(self, properties, initial_coordinates, initial_velocities):
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
        self.nodes = Node1d(nbr_nodes, initial_coordinates, initial_velocities,
                            section=properties.geometric.section)
        # ---------------------------------------------
        # Cells creation
        # ---------------------------------------------
        nbr_cells = nbr_nodes - 1
        self.cells = Element1d(nbr_cells, properties)
        # ---------------------------------------------
        # Topology creation
        # ---------------------------------------------
        self.__topologie = Topology1D(nbr_nodes, nbr_cells)
        # ---------------------------------------------
        # Ruptured cells vector
        # ---------------------------------------------
        self.__ruptured_cells = np.zeros(self.cells.number_of_cells, dtype=np.bool, order='C')

    def computeNodesMasses(self):
        """ Nodal mass computation """
        vecteur_nb_noeuds_par_element = np.zeros([self.cells.number_of_cells, ], dtype=np.int, order='C')
        vecteur_nb_noeuds_par_element[:] = 2
        self.nodes.calculer_masse_wilkins(self.__topologie, self.cells.masse, vecteur_nb_noeuds_par_element)

    def computeNewNodesVelocities(self, delta_t):
        """
        Computation of nodes velocities at t+dt

        :var delta_t: time step
        :type delta_t: float
        """
        self.nodes.calculer_nouvo_vitesse(delta_t)

    def computeNewNodesCoordinates(self, delta_t):
        """
        Computation of nodes coordinates at t+dt

        :var delta_t: time step
        :type delta_t: float
        """
        self.nodes.calculer_nouvo_coord(delta_t)

    def computeCellsSizes(self):
        """
        Computation of cells sizes at t
        """
        self.cells.computeSize(self.__topologie, self.nodes.xt)

    def computeNewCellsSizes(self, delta_t):
        """
        Computation of cells sizes at t+dt
        """
        self.cells.computeNewSize(self.__topologie, self.nodes.xtpdt, delta_t)

    def computeNewCellsDensities(self):
        """
        Computation of cells densities at t+dt
        """
        self.cells.computeNewDensity()

    def computeNewCellsPressures(self):
        """
        Computation of cells pressure at t+dt
        """
        self.cells.computeNewPressure(mask=~self.__ruptured_cells)

    def computeNewCellsPseudoViscosities(self, delta_t):
        """
        Computation of cells pseudoviscosities at t+dt

        :var delta_t: time step
        :type delta_t: float
        """
        self.cells.computeNewPseudo(delta_t)

    def computeNewNodesForces(self):
        """
        Computation of nodes forces at t+dt
        """
        self.nodes.calculer_nouvo_force(self.__topologie, self.cells.pressure.new_value, self.cells.pseudo.new_value)

    def increment(self):
        """
        Moving to next time step
        """
        self.nodes.incrementer()
        self.cells.incrementVariables()

    def computeNewTimeStep(self):
        """
        Computation of new time step
        """
        self.cells.computeNewTimeStep()
        return self.cells.dt.min()

    def applyPressure(self, surface, pressure):
        """
        Apply a given pressure on left or right boundary

        :var surface: name of the surface where pressure has to be imposed
        :var pressure: value of the pressure to impose
        :type surface: str ('left' | 'right')
        :type pressure: float
        """
        if surface.lower() not in ("left", "right"):
            raise(ValueError("One dimensional mesh : only 'left' or 'right' boundaries are possibles!"))
        if (surface.lower() == 'left'):
            self.nodes.applyPressure(0, pressure)
        else:
            self.nodes.applyPressure(-1, -pressure)

    def getRupturedCells(self, rupture_criterion):
        """
        Find the cells where the rupture criterion is checked and store them

        :var rupture_criterion: rupture criterion
        :type rupture_criterion: RuptureCriterion
        """
        self.__ruptured_cells = np.logical_or(self.__ruptured_cells, rupture_criterion.checkCriterion(self.cells))

    def applyRuptureTreatment(self, treatment):
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
        return self.cells.getCoordinates(self.cells.number_of_cells, self.__topologie, self.nodes.xt)

    @property
    def pressure_field(self):
        """
        Pressure field
        """
        return self.cells.pressure.current_value

    @property
    def density_field(self):
        """
        Density field
        """
        return self.cells.density.current_value

    @property
    def energy_field(self):
        """
        Internal energy field
        """
        return self.cells.energy.current_value

    @property
    def pseudoviscosity_field(self):
        """
        Pseudoviscosity field
        """
        return self.cells.pseudo.current_value


