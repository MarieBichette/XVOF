# -*- coding: utf-8 -*-
"""
Base class for one dimensional enriched mesh
"""

import numpy as np
from xfv.src.cell.one_dimension_enriched_cell_Hansbo import OneDimensionHansboEnrichedCell
from xfv.src.node.one_dimension_enriched_node_Hansbo import OneDimensionHansboEnrichedNode
from xfv.src.data.data_container import DataContainer
from xfv.src.mesh.topology1d import Topology1D
from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.utilities.profilingperso import timeit_file
from xfv.src.mass_matrix.one_dimension_mass_matrix import OneDimensionMassMatrix
from xfv.src.contact.contact import ContactModel

# noinspection PyArgumentList
class Mesh1dEnriched(object):  # pylint:disable=too-many-instance-attributes, too-many-public-methods
    """
    This class defines a one dimensional mesh with potential enrichment
    """
    # noinspection PyArgumentList
    def __init__(self, initial_coordinates, initial_velocities, enrichment_type):
        if np.shape(initial_coordinates) != np.shape(initial_velocities):
            message = "Initial velocity and coordinates vector doesn't have the same shape!"
            raise ValueError(message)
        if np.shape(initial_coordinates)[1] != 1:
            message = ("""A 1D mesh must have one dimensional vector which is not the case"""
                       """ for initial coordinates vector!""")
            raise ValueError(message)

        self.enrichment_type = enrichment_type

        # ---------------------------------------------
        # Nodes creation
        # ---------------------------------------------
        nbr_nodes = np.shape(initial_coordinates)[0]
        if self.enrichment_type == "Hansbo":
            self.nodes = OneDimensionHansboEnrichedNode(nbr_nodes, initial_coordinates,
                                                        initial_velocities,
                                                        section=DataContainer().geometric.section)
        else:
            self.nodes = OneDimensionHansboEnrichedNode(nbr_nodes, initial_coordinates,
                                                        initial_velocities,
                                                        section=DataContainer().geometric.section)

        # ---------------------------------------------
        # Cells creation
        # ---------------------------------------------
        nbr_cells = nbr_nodes - 1
        if self.enrichment_type == "Hansbo":
            self.cells = OneDimensionHansboEnrichedCell(nbr_cells)
        else:
            self.cells = OneDimensionHansboEnrichedCell(nbr_cells)

        # ----------------------------------------------
        # Mass Matrix creation
        # ----------------------------------------------
        self.mass_matrix = OneDimensionMassMatrix(nbr_nodes, correction_on_last_cells=None)
        # self.mass_matrix = OneDimensionMassMatrix(nbr_nodes, correction_on_last_cells="hansbo")

        # ---------------------------------------------
        # Topology creation
        # ---------------------------------------------
        self.__topology = Topology1D(nbr_nodes, nbr_cells)
        self.nb_nodes_per_cell = np.zeros([self.cells.number_of_cells, ], dtype=np.int, order='C')
        self.nb_nodes_per_cell[:] = 2

        # ---------------------------------------------
        # Ruptured cells vector
        # ---------------------------------------------
        self.__ruptured_cells = np.zeros(self.cells.number_of_cells, dtype=np.bool, order='C')
        self.__plastic_cells = np.zeros(self.cells.number_of_cells, dtype=np.bool, order='C')

        # Initialize cell fields
        self.cells.initialize_cell_fields(self.nodes.nodes_in_target,
                                          self.nodes.nodes_in_projectile,
                                          self.__topology)
        self.mask_last_nodes_of_ref = None

    @property
    def topology(self):
        """
        Retourne l'objet de topologie (pour qu'il soit accesible pour les tests unitaires)
        """
        return self.__topology

    def compute_cells_masses(self):
        """
        Cell mass computation
        """
        self.cells.compute_mass()

    def compute_nodes_masses(self):
        """
        Nodal mass computation
        """
        self.mass_matrix.compute_mass_matrix(self.__topology, self.cells.mass,
                                             self.nb_nodes_per_cell)

        if self.mass_matrix.correction_on_cell_500 is not None:
            print 'Matrix correction on last cells compatible with {} analyis'.format(
                self.mass_matrix.correction_on_cell_500)
            # identifier derniers éléments de la barre de référence
            self.mask_last_nodes_of_ref = np.zeros(
                [self.nodes.number_of_nodes], dtype=bool)
            self.mask_last_nodes_of_ref[-2] = True
            self.mask_last_nodes_of_ref[-1] = True
            self.mass_matrix.compute_correction_mass_matrix_for_cell_500(
                self.cells.mass, self.mask_last_nodes_of_ref, self.__topology)

    @timeit_file("/tmp/profil_xfv.src.txt")
    def compute_new_nodes_velocities(self, delta_t):
        """
        Computation of nodes velocities at t+dt
        :var delta_t: float, time step
        """
        # ddl classiques (loin de l'enrichissement)
        self.nodes.compute_new_velocity(
            delta_t, self.nodes.enrichment_not_concerned,
            self.mass_matrix.inverse_mass_matrix[self.nodes.enrichment_not_concerned])

        if self.mass_matrix.correction_on_cell_500 is not None:
            inv_mass_matrix_correction = self.mass_matrix.inverse_correction_mass_matrix
            # on applique la correction sur les derniers �l�ments de la matrice de r�f�rence :
            self.nodes.apply_correction_reference_bar(
                delta_t, inv_mass_matrix_correction,
                self.mass_matrix.inverse_mass_matrix[self.mask_last_nodes_of_ref],
                mask=self.mask_last_nodes_of_ref)
            #############
        for disc in Discontinuity.discontinuity_list():
            # Calcul des nouvelles matrices de masse enrichies pour les nouvelles discontinuit�s
            if not disc.mass_matrix_updated:
                # Construction de la matrice masse enrichie et de son inverse
                disc.mass_matrix_enriched.compute_enriched_mass_matrix(disc, self.__topology,
                                                                       self.cells.mass)
                if self.enrichment_type == "Hansbo":
                    disc.mass_matrix_enriched.assemble_enriched_mass_matrix(
                        "_enriched_mass_matrix_left_part", "_enriched_mass_matrix_right_part")
                    # réarrangement de la matrice de masse pour avoir
                    # strucure (classiq/enr/couplage)
                    disc.mass_matrix_enriched.rearrange_dof_in_inv_mass_matrix()
                disc.mass_matrix_enriched.print_enriched_mass_matrix()
                disc.has_mass_matrix_been_computed()

            # Calcul des vitesses ddl classique et enrichi � l'endroit de l'enrichissement
            self.nodes.compute_new_velocity(
                delta_t, disc.mask_disc_nodes,
                disc.mass_matrix_enriched.inverse_enriched_mass_matrix_classic_dof)
            self.nodes.compute_additional_dof_new_velocity(
                delta_t, disc.mass_matrix_enriched.inverse_enriched_mass_matrix_enriched_dof)

            # Couplage entre ddl classiques et enrichis
            self.nodes.coupled_enrichment_terms_compute_new_velocity(
                delta_t, disc.mass_matrix_enriched.inverse_enriched_mass_matrix_coupling_dof)

        self.nodes.compute_complete_velocity_field()

    def compute_contact(self, delta_t):
        """
        Compute the contact force to be applied to ensure non penetration of the
        discontinuities boundaries
        :param delta_t : time step, float
        """
        for disc in Discontinuity.discontinuity_list():
            contact = ContactModel(disc)
            contact.check_contact(self.nodes.xtpdt)
            if contact.has_contact:
                contact.compute_contact(self.nodes.upundemi, delta_t)
                contact.apply_contact(self.nodes.upundemi, delta_t)

    @timeit_file("/tmp/profil_xfv.src.txt")
    def compute_new_nodes_coordinates(self, delta_t):
        """
        Computation of nodes coordinates at t+dt

        :var delta_t: time step
        :type delta_t: float
        """
        self.nodes.compute_new_coodinates(delta_t)
        self.nodes.enriched_nodes_compute_new_coordinates(delta_t)

    @timeit_file("/tmp/profil_xfv.src.txt")
    def compute_cells_sizes(self):
        """
        Computation of cells sizes at t
        """
        self.cells.compute_size(self.__topology, self.nodes.xt)

    @timeit_file("/tmp/profil_xfv.src.txt")
    def compute_new_cells_sizes(self, delta_t):
        """
        Computation of cells sizes at t+dt
        """
        self.cells.compute_new_size(self.__topology, self.nodes.xtpdt, self.cells.classical)
        self.cells.compute_enriched_elements_new_part_size(delta_t, self.nodes.upundemi)
        self.nodes.compute_discontinuity_opening()

    @timeit_file("/tmp/profil_xfv.src.txt")
    def compute_new_cells_densities(self):
        """
        Computation of cells densities at t+dt
        """
        self.cells.compute_new_density(self.cells.classical)

        self.cells.compute_enriched_elements_new_density()

    @timeit_file("/tmp/profil_xfv.src.txt")
    def compute_new_cells_pressures(self, dt):  # dt name is ok pylint: disable=invalid-name
        """
        Computation of cells pressure at t+dt
        :var dt: time step
        :type dt: float
        """
        self.cells.compute_new_pressure(
            np.logical_and(self.cells.classical, ~self.__ruptured_cells), dt=dt)

        self.cells.compute_enriched_elements_new_pressure(dt)

    @timeit_file("/tmp/profil_xfv.src.txt")
    def compute_new_cells_pseudo_viscosity(self, delta_t):
        """
        Computation of cells artificial viscosity at t+dt
        :var delta_t: time step
        :type delta_t: float
        """
        self.cells.compute_new_pseudo(delta_t, self.cells.classical)
        self.cells.compute_enriched_elements_new_pseudo(delta_t)

    @timeit_file("/tmp/profil_xfv.src.txt")
    def compute_new_nodes_forces(self):
        """
        Computation of nodes forces at t+dt
        """
        self.nodes.compute_new_force(self.__topology, self.cells.stress_xx)
        self.nodes.compute_enriched_nodes_new_force(self.cells.stress_xx)

    @timeit_file("/tmp/profil_xfv.src.txt")
    def compute_new_cohesive_forces(self):
        """
        Computation of cohesive forces at t+dt
        """
        if DataContainer().material_target.failure_model.failure_treatment == "Enrichment":
            if DataContainer().material_target.damage_model.cohesive_model is not None:
                self.nodes.compute_enriched_nodes_cohesive_forces()

    @timeit_file("/tmp/profil_xfv.src.txt")
    def increment(self):
        """
        Moving to next time step
        """
        self.nodes.increment()
        self.cells.increment_variables()
        for disc in Discontinuity.discontinuity_list():
            disc.additional_dof_increment()

    @timeit_file("/tmp/profil_xfv.src.txt")
    def compute_deviator_elasticity(self, delta_t, mask):
        """
        Compute the deviatoric part of stress tensor
        :param delta_t : float, time step staggered
        :param mask: array of bool to select cells of interest
        """
        # Sert à identifier si on est dans le  projectile ou dans la cible
        mask = np.logical_and(
            mask, self.cells.classical)  # pylint: disable=assignment-from-no-return
        self.cells.compute_deviatoric_stress_tensor(mask, self.__topology,
                                                    self.nodes.xtpdt, self.nodes.upundemi, delta_t)
        self.cells.compute_enriched_deviatoric_stress_tensor(self.nodes.xtpdt,
                                                             self.nodes.upundemi, delta_t)

    @timeit_file("/tmp/profil_xfv.src.txt")
    def assemble_complete_stress_tensor(self):
        """
        Assembling pressure and stress deviator
        """
        self.cells.compute_complete_stress_tensor(self.cells.classical)
        self.cells.compute_enriched_stress_tensor()

    @timeit_file("/tmp/profil_xfv.src.txt")
    def compute_new_time_step(self):
        """
        Computation of new time step
        """
        # self.cells.compute_new_time_step(self.cells.classical)
        # self.cells.compute_enriched_elements_new_time_step()
        # return self.cells.dt.min()

        initial_time_step = DataContainer().time.initial_time_step
        dt = initial_time_step  # dt name is ok pylint: disable=invalid-name
        reduction_factor = DataContainer().time.time_step_reduction_factor_for_failure
        if reduction_factor is not None:
            if self.cells.enriched.any():
                dt = dt/reduction_factor  # dt name is ok pylint: disable=invalid-name
        return dt

    @timeit_file("/tmp/profil_xfv.src.txt")
    def apply_pressure(self, surface, pressure):
        """
        Apply a given pressure on left or right boundary
        :var surface: name of the surface where pressure has to be imposed
        :var pressure: value of the pressure to impose
        :type surface: str ('left' | 'right')
        :type pressure: float
        """
        if surface.lower() not in ("left", "right"):
            raise(ValueError("""One dimensional mesh : only """
                             """'left' or 'right' boundaries are possibles!"""))
        if surface.lower() == 'left':
            self.nodes.apply_pressure(0, pressure)
        else:
            self.nodes.apply_pressure(-1, -pressure)

    @timeit_file("/tmp/profil_xfv.src.txt")
    def apply_velocity_boundary_condition(self, surface, velocity):
        """
        Apply a given velocity on left or right boundary
        """
        if surface.lower() not in ("left", "right"):
            raise ValueError("""One dimensional mesh : only """
                             """'left' or 'right' boundaries are possibles!""")
        if surface.lower() == 'left':
            self.nodes.apply_velocity_boundary_coundition(0, velocity)
        else:
            self.nodes.apply_velocity_boundary_coundition(-1, velocity)

    @timeit_file("/tmp/profil_xfv.src.txt")
    def get_ruptured_cells(self, rupture_criterion):
        """
        Find the cells where the rupture criterion is checked and store them
        :var rupture_criterion: rupture criterion
        :type rupture_criterion: RuptureCriterion
        """
        new_cracked_cells_in_target = rupture_criterion.check_criterion(self.cells)
        # correction car le projectile ne peut pas rompre
        new_cracked_cells_in_target[self.cells.cell_in_projectile] = False
        self.__ruptured_cells = \
            np.logical_or(self.__ruptured_cells,
                          new_cracked_cells_in_target) # pylint: disable=assignment-from-no-return

    @timeit_file("/tmp/profil_xfv.src.txt")
    def get_plastic_cells(self, plastic_criterion, mask):
        """
        Find the cells where the plasticity criterion is checked and store them
        :var plastic_criterion: plastic criterion
        :type plastic_criterion: PlasticityCriterion
        :param mask: array of bool to select cells of interest
        """
        self.__plastic_cells[mask] = plastic_criterion.check_criterion(self.cells)[mask]

        for disc in Discontinuity.discontinuity_list():
            disc.plastic_cells = plastic_criterion.check_criterion_on_right_part_cells(disc)


    @timeit_file("/tmp/profil_xfv.src.txt")
    def apply_rupture_treatment(self, treatment, time):
        """
        Apply the rupture treatment on the cells enforcing the rupture criterion
        :var treatment: rupture treatment
        :type treatment: RuptureTreatment
        """
        treatment.apply_treatment(self.cells, self.__ruptured_cells,
                                  self.nodes, self.__topology, time)

    @timeit_file("/tmp/profil_xfv.src.txt")
    def apply_plasticity_treatment(self, dt):  # dt is ok. pylint: disable=invalid-name
        """
        Apply plasticity treatment if criterion is activated :
        - compute plastic deviatoric stress tensor
        - compute plastic strain rate
        - compute yield stress
        :param dt : float, time step
        """
        # La méthode apply_plastic_corrector_on_deviatoric_stress_tensor modifie
        # la variable dev_stress_new et doit donc être appelée à la fin de
        # l'étape du calcul de plasticité pour conserver la prédiction élastique dans
        # le calcul du taux de déformation plastique, plasticité cumulée, ...

        # Cells plastic classical
        mask = np.logical_and(
            self.cells.classical, self.__plastic_cells) # pylint: disable=assignment-from-no-return
        self.cells.compute_yield_stress()
        self.cells.compute_plastic_strain_rate_tensor(mask, dt)
        self.cells.compute_equivalent_plastic_strain_rate(mask, dt)
        self.cells.apply_plastic_corrector_on_deviatoric_stress_tensor(mask)

        # Cells plastic enriched
        self.cells.compute_enriched_yield_stress()
        self.cells.compute_enriched_plastic_strain_rate(self.__plastic_cells, dt)
        self.cells.compute_enriched_equivalent_plastic_strain_rate(self.__plastic_cells, dt)
        self.cells.apply_plastic_correction_on_enriched_deviatoric_stress_tensor(
            self.__plastic_cells)

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
        return self.nodes.xt

    @property
    def cells_coordinates(self):
        """
        Cells coordinates (coordinates of cells centers)
        """
        # Pour reconstruire le champ de coordonnées des cells,
        # les ruptured cells des discontinuités doivent être
        # triées par cell id pour savoir comment gérer le décalage
        modified_coord = np.zeros([len(Discontinuity.discontinuity_list()), 3])
        # modified_coord est un array qui contient ruptured_cell_id, left_size, right_size
        for disc in Discontinuity.discontinuity_list():
            modified_coord[disc.label - 1, 0] = int(disc.ruptured_cell_id)
            modified_coord[disc.label - 1, 1] = \
                self.nodes.xt[disc.mask_in_nodes] + disc.left_part_size.current_value / 2.
            modified_coord[disc.label - 1, 2] = \
                self.nodes.xt[disc.mask_out_nodes] - disc.right_part_size.current_value / 2.
        modified_coord = np.sort(modified_coord, 0)

        res = self.cells.get_coordinates(self.cells.number_of_cells, self.__topology, self.nodes.xt)
        # On sépare  les deux étapes de construction de cell_coordinates
        # pour ne pas écraser les résultats au fur et à mesure
        # Etape 1 : modification des longueurs de cell qui sont rompues :
        # taille non rompue ->taille gauche
        ligne_indice_cell = 0
        for indice_cell_rompue in modified_coord[:, 0]:
            res[int(indice_cell_rompue)] = modified_coord[ligne_indice_cell, 1]
            ligne_indice_cell += 1
        # Etape 2 : insertion des tailles droites
        ligne_indice_cell = 0
        for indice_cell_rompue in modified_coord[:, 0]:
            res = np.insert(res, int(indice_cell_rompue) + ligne_indice_cell + 1,
                            modified_coord[ligne_indice_cell, 2])
            ligne_indice_cell += 1
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

    @property
    def deviatoric_stress_field(self):
        """
        Deviatoric stress field
        """
        return self.cells.deviatoric_stress_field

    @property
    def stress_xx_field(self):
        """
        First component of the Cauchy stress tensor
        """
        return self.cells.stress_xx_field
