# -*- coding: iso-8859-1 -*-
"""
Base class for one dimensional mesh
"""
import os

import numpy as np

from xvof.cell.one_dimension_enriched_cell import OneDimensionEnrichedCell
from xvof.data.data_container import DataContainer
from xvof.mass_matrix.one_dimension_mass_matrix import OneDimensionMassMatrix
from xvof.mass_matrix.one_dimension_enriched_mass_matrix import OneDimensionEnrichedMassMatrix
from xvof.mesh.topology1d import Topology1D
from xvof.node.one_dimension_enriched_node import OneDimensionEnrichedNode
from xvof.discontinuity.discontinuity import Discontinuity
from xvof.utilities.profilingperso import timeit_file


class Mesh1dEnriched(object):
    """
    This class defines a one dimensional mesh with potential enrichment
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
        self.nodes = OneDimensionEnrichedNode(nbr_nodes, initial_coordinates, initial_velocities,
                                              section=DataContainer().geometric.section)
        # ---------------------------------------------
        # Cells creation
        # ---------------------------------------------
        nbr_cells = nbr_nodes - 1
        self.cells = OneDimensionEnrichedCell(nbr_cells)
        self.nb_nodes_per_cell = np.zeros([self.cells.number_of_cells, ], dtype=np.int, order='C')
        self.nb_nodes_per_cell[:] = 2
        # ---------------------------------------------
        # Topology creation
        # ---------------------------------------------
        self.__topology = Topology1D(nbr_nodes, nbr_cells)
        # ---------------------------------------------
        # Ruptured cells vector
        # ---------------------------------------------
        self.__ruptured_cells = np.zeros(self.cells.number_of_cells, dtype=np.bool, order='C')
        # ----------------------------------------------
        # Mass Matrix creation
        # ----------------------------------------------
        self.mass_matrix = OneDimensionMassMatrix(nbr_nodes, correction_3x3_on_cell_500=True)
        self.mass_matrix_enriched = OneDimensionEnrichedMassMatrix(lumped_matrix_classic_dof=False,
                                                                   lumped_matrix_enr_dof=False,
                                                                   analytical_inverse=False)

    def compute_cells_masses(self):
        """ Cell mass computation """
        self.cells.compute_mass()

    def compute_nodes_masses(self):
        """ node mass computation """
        self.mass_matrix.compute_mass_matrix(self.__topology, self.cells.mass, self.nb_nodes_per_cell)
        # import  ipdb ;ipdb.set_trace()
        if self.mass_matrix.correction_3x3_on_cell_500:
            self.mass_matrix.compute_3x3_mass_matrix_for_cell_500(self.cells.mass)
            self.mask_last_cells_of_ref = np.empty([self.nodes.number_of_nodes, 1], dtype=bool)
            self.mask_last_cells_of_ref[:] = False # mask pour identifier les deux derniers éléments de la barre de référence
            self.mask_last_cells_of_ref[-3] = True
            self.mask_last_cells_of_ref[-1] = True
            self.mask_last_cells_of_ref[-2] = True
            self.inverse_correction = self.mass_matrix.inverse_3x3_mass_matrix

    @timeit_file("/tmp/profil_xvof.txt")
    def compute_new_nodes_velocities(self, delta_t):
        """
        Computation of nodes velocities at t+dt

        :var delta_t: time step
        :type delta_t: float
        """
        # ddl classiques (loin de l'enrichissement )
        self.nodes.compute_new_velocity(delta_t, self.nodes.enrichment_not_concerned,
                                        self.mass_matrix.inverse_mass_matrix[self.nodes.enrichment_not_concerned])

        # on applique la correction sur le dernier élément de la matrice de référence :
        if self.mass_matrix.correction_3x3_on_cell_500:
            self.nodes.apply_correction_for_complete_mass_matrix_cell_500_ref(delta_t,
                                                                         self.inverse_correction,
                                                                         self.mass_matrix.inverse_mass_matrix
                                                                           [self.mask_last_cells_of_ref],
                                                                         mask=self.mask_last_cells_of_ref)

        if self.nodes.enriched.any():
            for disc in [d for d in Discontinuity.discontinuity_list() if not d.mass_matrix_updated]:
                # Construction de la matrice masse enrichie et de son inverse
                self.mass_matrix_enriched.compute_enriched_mass_matrix(self.__topology, self.cells.mass)
                self.mass_matrix_enriched.assemble_enriched_mass_matrix("__matrix_classic_dof", "__matrix_enr_dof",
                                                                        "__matrix_coupling")
                # self.mass_matrix_enriched.assemble_enriched_mass_matrix("__matrix_classic_dof", "__matrix_enr_dof")
                self.mass_matrix_enriched.print_enriched_mass_matrix()
                disc.hasMassMatrixBeenComputed()
            # Calcule des vitesses ddl classique et enrichi à l'endroit de l'enrichissement
            self.nodes.compute_new_velocity(delta_t, self.nodes.enrichment_concerned,
                                            self.mass_matrix_enriched.inverse_enriched_mass_matrix_classic_dof)
            self.nodes.enriched_nodes_compute_new_velocity(delta_t, self.nodes.enriched,
                                            self.mass_matrix_enriched.inverse_enriched_mass_matrix_enriched_dof)
            # Couplage entre ddl classiques et enrichis
            self.nodes.coupled_enrichment_terms_compute_new_velocity(delta_t,
                                            self.mass_matrix_enriched.inverse_enriched_mass_matrix_coupling_dof)
        self.nodes.compute_complete_velocity_field()

    @timeit_file("/tmp/profil_xvof.txt")
    def compute_new_nodes_coordinates(self, delta_t):
        """
        Computation of nodes coordinates at t+dt

        :var delta_t: time step
        :type delta_t: float
        """
        self.nodes.compute_new_coodinates(delta_t)
        self.nodes.enriched_nodes_compute_new_coordinates(delta_t)

    @timeit_file("/tmp/profil_xvof.txt")
    def compute_cells_sizes(self):
        """
        Computation of cells sizes at t
        """
        self.cells.compute_size(self.__topology, self.nodes.xt)


    @timeit_file("/tmp/profil_xvof.txt")
    def compute_new_cells_sizes(self, delta_t):
        """
        Computation of cells sizes at t+dt
        """
        self.cells.compute_new_size(self.__topology, self.nodes.xtpdt, mask=self.cells.classical)
        self.cells.compute_enriched_elements_new_part_size(delta_t, self.__topology, self.nodes.upundemi_enriched,
                                                           self.nodes.upundemi)

    @timeit_file("/tmp/profil_xvof.txt")
    def compute_new_cells_densities(self):
        """
        Computation of cells densities at t+dt
        """
        self.cells.compute_new_density(self.cells.classical)
        self.cells.compute_enriched_elements_new_density()

    @timeit_file("/tmp/profil_xvof.txt")
    def compute_new_cells_pressures(self):
        """
        Computation of cells pressure at t+dt
        """
        self.cells.compute_new_pressure(mask=np.logical_and(self.cells.classical, ~self.__ruptured_cells))
        self.cells.compute_enriched_elements_new_pressure()

    @timeit_file("/tmp/profil_xvof.txt")
    def compute_new_cells_pseudo_viscosity(self, delta_t):
        """
        Computation of cells artificial viscosity at t+dt

        :var delta_t: time step
        :type delta_t: float
        """
        self.cells.compute_new_pseudo(delta_t, mask=self.cells.classical)
        self.cells.compute_enriched_elements_new_pseudo(delta_t)

    @timeit_file("/tmp/profil_xvof.txt")
    def compute_new_nodes_forces(self):
        """
        Computation of nodes forces at t+dt
        """
        self.nodes.compute_new_force(self.__topology, self.cells.pressure.new_value, self.cells.pseudo.new_value)
        self.nodes.enriched_nodes_compute_new_force(self.__topology, self.cells.pressure.new_value,
                                                    self.cells.pressure.new_enr_value, self.cells.pseudo.new_value,
                                                    self.cells.pseudo.new_enr_value)

    @timeit_file("/tmp/profil_xvof.txt")
    def increment(self):
        """
        Moving to next time step
        """
        self.nodes.increment()
        self.nodes.enriched_nodes_increment()
        self.cells.increment_variables()

    @timeit_file("/tmp/profil_xvof.txt")
    def compute_new_time_step(self):
        """
        Computation of new time step
        """
        self.cells.compute_new_time_step(self.cells.classical)
        self.cells.compute_enriched_elements_new_time_step()
        return self.cells.dt.min()

    @timeit_file("/tmp/profil_xvof.txt")
    def apply_pressure(self, surface, pressure):
        """
        Apply a given pressure on left or right boundary

        :var surface: name of the surface where pressure has to be imposed
        :var pressure: value of the pressure to impose
        :type surface: str ('left' | 'right')
        :type pressure: float
        """
        if surface.lower() not in ("left", "right"):
            raise(ValueError("One dimensional mesh : only 'left' or 'right' boundaries are possibles!"))
        if surface.lower() == 'left':
            self.nodes.apply_pressure(0, pressure)
        else:
            self.nodes.apply_pressure(-1, -pressure)

    @timeit_file("/tmp/profil_xvof.txt")
    def get_ruptured_cells(self, rupture_criterion):
        """
        Find the cells where the rupture criterion is checked and store them

        :var rupture_criterion: rupture criterion
        :type rupture_criterion: RuptureCriterion
        """
        self.__ruptured_cells = np.logical_or(self.__ruptured_cells, rupture_criterion.checkCriterion(self.cells))

    @timeit_file("/tmp/profil_xvof.txt")
    def apply_rupture_treatment(self, treatment):
        """
        Apply the rupture treatment on the cells enforcing the rupture criterion

        :var treatment: rupture treatment
        :type treatment: RuptureTreatment
        """
        treatment.applyTreatment(self.cells, self.__ruptured_cells, self.nodes, self.__topology)

    @property
    def velocity_field(self):
        """
        Node velocity field
        """
        return self.nodes.velocity_field

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
        res = self.cells.get_coordinates(self.cells.number_of_cells, self.__topology, self.nodes.xt)
        for i in xrange(self.cells.number_of_cells):
            if self.cells.enriched[i]:
                nodes_index = self.__topology.getNodesBelongingToCell(i)
                res[i] = self.nodes.xtpdt[nodes_index][0] + self.cells.left_size.new_value[i] / 2.
                # import ipdb ; ipdb.set_trace()
                res = np.insert(res, i + 1, self.nodes.xtpdt[nodes_index][1] - self.cells.right_size.new_value[i] / 2.,
                                axis=0)
        return res 

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
