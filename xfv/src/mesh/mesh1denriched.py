# -*- coding: utf-8 -*-
"""
Base class for one dimensional enriched mesh
"""

import numpy as np
from xfv.src.cell.one_dimension_enriched_cell_Hansbo import OneDimensionHansboEnrichedCell
from xfv.src.node.one_dimension_enriched_node_Hansbo import OneDimensionHansboEnrichedNode
from xfv.src.data.data_container import DataContainer
from xfv.src.data.enriched_mass_matrix_props import ConsistentMassMatrixProps
from xfv.src.mesh.topology1d import Topology1D
from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.mass_matrix.one_dimension_mass_matrix import OneDimensionMassMatrix
from xfv.src.contact.contact_base import ContactBase


# noinspection PyArgumentList
class Mesh1dEnriched:  # pylint:disable=too-many-instance-attributes, too-many-public-methods
    """
    This class defines a one dimensional mesh with potential enrichment
    """
    # noinspection PyArgumentList
    def __init__(self, initial_coordinates, initial_velocities):
        """
        Construction of the mesh
        :param initial_coordinates:
        :param initial_velocities:
        """
        self.data = DataContainer()  # pylint: disable=no-value-for-parameter
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
        self.nodes = OneDimensionHansboEnrichedNode(nbr_nodes, initial_coordinates,
                                                    initial_velocities,
                                                    section=self.data.geometric.section)

        # ---------------------------------------------
        # Cells creation
        # ---------------------------------------------
        nbr_cells = nbr_nodes - 1
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

        # ---------------------------------------------
        # Cohesive zone model initialisation
        # ---------------------------------------------
        self.cohesive_zone_model = None
        target_dmg_model = self.data.material_target.damage_model
        if target_dmg_model is not None:
            self.cohesive_zone_model = target_dmg_model.cohesive_model.build_cohesive_model_obj()
        if (self.data.material_target.failure_model.failure_treatment != "Enrichment") and \
                (self.cohesive_zone_model is not None):
            print("No cohesive model is allowed if failure treatment is not Enrichment")
            self.cohesive_zone_model = None

        # ---------------------------------------------
        # Contact model initialisation (contact between discontinuities boundaries)
        # ---------------------------------------------
        if self.data.material_target.contact_model is not None:
            self.contact_model: ContactBase = \
                self.data.material_target.contact_model.contact_model.build_contact_obj()
        else:
            self.contact_model = None

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
            print('Matrix correction on last cells compatible with {} analyis'.format(
                self.mass_matrix.correction_on_cell_500))
            # Identify the last elements of the reference bar
            self.mask_last_nodes_of_ref = np.zeros(
                [self.nodes.number_of_nodes], dtype=bool)
            self.mask_last_nodes_of_ref[-2] = True
            self.mask_last_nodes_of_ref[-1] = True
            self.mass_matrix.compute_correction_mass_matrix_for_cell_500(
                self.cells.mass, self.mask_last_nodes_of_ref, self.__topology)

    def compute_new_nodes_velocities(self, delta_t: float):
        """
        Computation of nodes velocities at t+dt
        :var delta_t: time step
        """
        self._compute_velocities_for_enrichment_not_concerned_nodes(delta_t)
        # Compute velocity for enriched nodes
        for disc in Discontinuity.discontinuity_list():
            # Compute mass matrix for newly created discontinuities
            if not disc.mass_matrix_updated:
                self._compute_discontinuity_mass_matrix(disc)
            # Compute classical and enriched ddl velocity of enriched nodes
            self._compute_velocities_for_disc(disc, delta_t)
        self.nodes.compute_complete_velocity_field()

    def _compute_velocities_for_enrichment_not_concerned_nodes(self, delta_t: float):
        """
        Compute classical ddl (far from enrichment = enrichment not concerned)
        :param delta_t: time step
        """
        # Compute classical ddl (far from enrichment = enrichment not concerned)
        self.nodes.compute_new_velocity(
            delta_t, self.nodes.enrichment_not_concerned,
            self.mass_matrix.inverse_mass_matrix[self.nodes.enrichment_not_concerned])

        if self.mass_matrix.correction_on_cell_500 is not None:
            # Apply some correction to mimic a consistent mass matrix on the last cells of
            # the reference bar
            inv_mass_matrix_correction = self.mass_matrix.inverse_correction_mass_matrix
            self.nodes.apply_correction_reference_bar(
                delta_t, inv_mass_matrix_correction,
                self.mass_matrix.inverse_mass_matrix[self.mask_last_nodes_of_ref],
                mask=self.mask_last_nodes_of_ref)

    def _compute_velocities_for_disc(self, disc: Discontinuity, delta_t: float):
        """
        Compute the new node velocities for the nodes belonging to a given discontinuity
        :param disc :current discontinuity
        :param delta_t: time step
        """
        # Compute classical ddl velocity of enriched nodes
        self.nodes.compute_new_velocity(
            delta_t, disc.mask_disc_nodes,
            disc.mass_matrix_enriched.inverse_enriched_mass_matrix_classic_dof)
        # Compute enriched ddl velocity of enriched nodes
        self.nodes.compute_additional_dof_new_velocity(
            delta_t, disc.mass_matrix_enriched.inverse_enriched_mass_matrix_enriched_dof)

        if type(self.data.material_target.failure_model.lump_mass_matrix) == \
                ConsistentMassMatrixProps:
            # Compute the contribution of classical ddl on enriched ddl and the reverse
            # (out of the diagonal terms of the mass matrix)
            self.nodes.coupled_enrichment_terms_compute_new_velocity(
                delta_t, disc.mass_matrix_enriched.inverse_enriched_mass_matrix_coupling_dof)

    def _compute_discontinuity_mass_matrix(self, disc: Discontinuity):
        """
        Compute the mass matrix of a newly created discontinuity
        :param disc: Discontinuity
        """
        disc.mass_matrix_enriched.compute_enriched_mass_matrix(disc, self.__topology,
                                                               self.cells.mass)
        disc.mass_matrix_enriched.assemble_enriched_mass_matrix(
            "_enriched_mass_matrix_left_part", "_enriched_mass_matrix_right_part")
        # Arrange mass matrix to get a structure : classical / enriched dof
        disc.mass_matrix_enriched.rearrange_dof_in_inv_mass_matrix()
        disc.has_mass_matrix_been_computed()

    def apply_contact_correction(self, delta_t: float):
        """
        Compute the contact force to be applied to ensure non penetration of the
        discontinuities boundaries
        :param delta_t : time step
        """
        if self.contact_model is not None:
            # Theoretically, we should consider a global resolution of contact in all
            # discontinuities. Here they are treated one after another. Better than nothing but
            # may cause instabilities
            for disc in Discontinuity.discontinuity_list():
                contact_force = self.contact_model.compute_contact_force(self.nodes.upundemi,
                                                                         disc, delta_t)
                if contact_force != 0.:
                    # Divide the contact "force" on the nodal forces
                    self.nodes.apply_force_on_discontinuity_boundaries(disc, contact_force)

                    # Reinitialize the kinematics that lead to contact
                    self.nodes.reinitialize_kinematics_after_contact(disc)
                    disc.reinitialize_kinematics_after_contact()

                    # Apply correction on the velocity field (only on disc nodes)
                    self._compute_velocities_for_disc(disc, delta_t)

                    # Theoretically, we should apply the velocity boundary condition here,
                    # but it is really not convenient to do this and fracture is not supposed
                    # to occur on the boundary cells. Thus, no boundary conditions is applied

                    # Apply correction on the node coordinates (only on disc nodes)
                    self.nodes.compute_new_coodinates(disc.mask_disc_nodes, delta_t)  # classical
                    self.nodes.enriched_nodes_compute_new_coordinates(disc, delta_t)  # enriched
                    # Update discontinuity opening
                    disc.compute_discontinuity_new_opening(self.nodes.xtpdt)  # opening

    def compute_new_nodes_coordinates(self, delta_t: float):
        """
        Computation of nodes coordinates at t+dt
        :param delta_t: time step
        """
        mask_all_nodes = np.ones([self.nodes.number_of_nodes], dtype=bool)
        self.nodes.compute_new_coodinates(mask_all_nodes, delta_t)
        for disc in Discontinuity.discontinuity_list():
            self.nodes.enriched_nodes_compute_new_coordinates(disc, delta_t)
            # Update discontinuity opening
            disc.compute_discontinuity_new_opening(self.nodes.xtpdt)

    def compute_cells_sizes(self):
        """
        Computation of cells sizes at t
        """
        self.cells.compute_size(self.__topology, self.nodes.xt)

    def compute_new_cells_sizes(self, delta_t):
        """
        Computation of cells sizes at t+dt
        """
        self.cells.compute_new_size(self.__topology, self.nodes.xtpdt, self.cells.classical)
        self.cells.compute_enriched_elements_new_part_size(delta_t, self.nodes.upundemi)

    def compute_new_cells_densities(self):
        """
        Computation of cells densities at t+dt
        """
        self.cells.compute_new_density(self.cells.classical)
        self.cells.compute_enriched_elements_new_density()

    def compute_new_cells_pressures(self, delta_t: float):
        """
        Computation of cells pressure at t+dt
        :var delta_t: time step
        """
        # all classical cells + left part of enriched cells :
        self.cells.compute_new_pressure(~self.__ruptured_cells, dt=delta_t)
        # right part of enriched cells
        if self.cells.enriched.any():
            # test if any enriched cell in order to avoid error in the Newton initialization
            self.cells.compute_enriched_elements_new_pressure(delta_t)

    def compute_new_cells_pseudo_viscosity(self, delta_t):
        """
        Computation of cells artificial viscosity at t+dt
        :var delta_t: time step
        :type delta_t: float
        """
        self.cells.compute_new_pseudo(delta_t, self.cells.classical)
        self.cells.compute_enriched_elements_new_pseudo(delta_t)

    def compute_new_nodes_forces(self):
        """
        Computation of nodes forces at t+dt
        """
        self.nodes.compute_new_force(self.__topology, self.cells.stress_xx, self.cells.classical)
        self.nodes.compute_enriched_nodes_new_force(self.cells.stress_xx,
                                                    self.cells.additional_dof_stress_xx)

    def compute_new_cohesive_forces(self):
        """
        Computation of cohesive forces at t+dt
        """
        if self.cohesive_zone_model is not None:
            self.nodes.compute_enriched_nodes_cohesive_forces(self.cohesive_zone_model)

    def increment(self):
        """
        Moving to next time step
        """
        self.nodes.increment()
        self.cells.increment_variables()
        self.cells.cell_additional_dof_increment()  # enriched cell variables
        for disc in Discontinuity.discontinuity_list():
            disc.additional_dof_increment()  # enriched node variables

    def apply_elasticity(self, delta_t, shear_modulus_model, mask_material):
        """
        Compute the deviatoric part of stress tensor
        :param delta_t : float, time step staggered
        :param shear_modulus_model: model to compute the shear modulus
        :param mask_material: array of bool to select cells of interest
        """
        # Sert à identifier si on est dans le  projectile ou dans la cible
        mask = np.logical_and(
            mask_material, self.cells.classical)  # pylint: disable=assignment-from-no-return

        # Update the shear modulus (same calculation classic and enr left => mask_material)
        self.cells.compute_shear_modulus(shear_modulus_model, mask_material)
        if mask_material.any():
            # avoid calculation if mask_material = projectile because no enrichment in projectile
            self.cells.compute_enriched_shear_modulus(shear_modulus_model)

        # Compute the deviatoric stress tensor (only classical cells)
        self.cells.compute_deviatoric_stress_tensor(mask, self.__topology,
                                                    self.nodes.xtpdt, self.nodes.upundemi, delta_t)
        # Compute dev stress tensor for enr cells (left and right parts)
        # (full enr ddl because coord and velocities are required for strain rate computation)
        if mask_material.any():
            # avoid calculation if mask_material = projectile because no enrichment in projectile
            self.cells.compute_enriched_deviatoric_stress_tensor(self.nodes.xtpdt,
                                                                 self.nodes.upundemi, delta_t)

    def assemble_complete_stress_tensor(self):
        """
        Assembling pressure and stress deviator
        """
        self.cells.compute_complete_stress_tensor()
        if self.cells.enriched.any():
            self.cells.compute_enriched_stress_tensor()

    def compute_new_time_step(self):
        """
        Computation of new time step
        """
        if not self.data.time.is_time_step_constant:
            self.cells.compute_new_time_step(self.cells.classical)
            self.cells.compute_enriched_elements_new_time_step()
            dt = self.cells.dt.min()  # dt name is ok pylint: disable=invalid-name
        else:
            initial_time_step = self.data.time.initial_time_step
            dt = initial_time_step  # dt name is ok pylint: disable=invalid-name

        reduction_factor = self.data.time.time_step_reduction_factor_for_failure
        if reduction_factor is not None:
            if self.cells.enriched.any():
                dt = dt/reduction_factor  # dt name is ok pylint: disable=invalid-name
        return dt

    def apply_pressure(self, surface, pressure):
        """
        Apply a given pressure on left or right boundary
        :var surface: name of the surface where pressure has to be imposed
        :var pressure: value of the pressure to impose
        :type surface: str ('left' | 'right')
        :type pressure: float
        """
        if surface.lower() == 'left':
            self.nodes.apply_pressure(0, pressure)
        elif surface.lower() == 'right':
            self.nodes.apply_pressure(-1, -pressure)
        else:
            msg = "One dimensional mesh : only 'left' or 'right' boundaries are possibles!"
            raise ValueError(msg)

    def apply_velocity_boundary_condition(self, surface, velocity):
        """
        Apply a given velocity on left or right boundary
        """
        if surface.lower() == 'left':
            self.nodes.apply_velocity_boundary_coundition(0, velocity)
        elif surface.lower() == 'right':
            self.nodes.apply_velocity_boundary_coundition(-1, velocity)
        else:
            msg = "One dimensional mesh : only 'left' or 'right' boundaries are possibles!"
            raise ValueError(msg)

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
                          new_cracked_cells_in_target)  # pylint: disable=assignment-from-no-return

    def _get_plastic_cells(self, plastic_criterion, mask):
        """
        Find the cells where the plasticity criterion is checked and store them
        :var plastic_criterion: plastic criterion
        :type plastic_criterion: PlasticityCriterion
        :param mask: array of bool to select cells of interest
        """
        self.__plastic_cells[mask] = plastic_criterion.check_criterion(self.cells)[mask]
        self.cells.plastic_enr_cells[mask] = \
            plastic_criterion.check_criterion_on_right_part_cells(self.cells)[mask]

    def apply_rupture_treatment(self, treatment, time: float):
        """
        Apply the rupture treatment on the cells enforcing the rupture criterion
        :var treatment: rupture treatment
        :type treatment: RuptureTreatment
        :param time : simulation time
        """
        treatment.apply_treatment(self.cells, self.__ruptured_cells,
                                  self.nodes, self.__topology, time)

    def apply_plasticity(self, delta_t: float, yield_stress_model, plasticity_criterion,
                         mask_mesh: np.array):
        """
        Apply plasticity treatment if criterion is activated :
        - compute yield stress
        - tests plasticity criterion
        - compute plastic strain rate for plastic cells
        :param delta_t : time step
        :param yield_stress_model: model to compute the yield stress
        :param plasticity_criterion: model for the plasticity criterion
        :param mask_mesh : mask cells in projectile or target
        """
        # La méthode apply_plastic_corrector_on_deviatoric_stress_tensor modifie
        # la variable dev_stress_new et doit donc être appelée à la fin de
        # l'étape du calcul de plasticité pour conserver la prédiction élastique dans
        # le calcul du taux de déformation plastique, plasticité cumulée, ...

        # 1) Compute yield stress
        self.cells.compute_yield_stress(yield_stress_model, mask_mesh)
        if mask_mesh.any():
            self.cells.compute_enriched_yield_stress(yield_stress_model)

        # 2) Get plastic cells (verification of the plasticity criterion)
        # Criterion is tested for both classical and enriched cells
        self._get_plastic_cells(plasticity_criterion, mask_mesh)

        # Get plastic cells either in projectile or in target
        mask = np.logical_and(mask_mesh,
                              self.__plastic_cells)  # pylint: disable=assignment-from-no-return
        # 3) Plasticity treatment for classical plastic cells and left part of enriched cells
        self.cells.compute_plastic_strain_rate_tensor(mask, delta_t)
        self.cells.compute_equivalent_plastic_strain_rate(mask, delta_t)
        self.cells.apply_plastic_corrector_on_deviatoric_stress_tensor(mask)

        # 4) Plasticity treatment for enriched plastic cells (right part)
        self.cells.compute_enriched_plastic_strain_rate(mask_mesh, delta_t)
        self.cells.compute_enriched_equivalent_plastic_strain_rate(mask_mesh, delta_t)
        self.cells.apply_plastic_correction_on_enriched_deviatoric_stress_tensor(mask_mesh)

    @property
    def velocity_field(self) -> np.array:
        """
        Node velocity field
        """
        return self.nodes.velocity_field

    @property
    def nodes_coordinates(self) -> np.array:
        """
        Nodes coordinates
        """
        return self.nodes.xt

    @property
    def cells_coordinates(self) -> np.array:
        """
        Cells coordinates (coordinates of cells centers)
        """
        # Pour reconstruire le champ de coordonnées des cells,
        # les ruptured cells des discontinuités doivent être
        # triées par cell id pour savoir comment gérer le décalage
        modified_coord = np.zeros([len(Discontinuity.discontinuity_list()), 3])
        # modified_coord est un array qui contient ruptured_cell_id, left_size, right_size
        for disc in Discontinuity.discontinuity_list():
            enr_cell = int(disc.ruptured_cell_id)
            index = disc.label - 1
            modified_coord[index, 0] = enr_cell
            modified_coord[index, 1] = self.nodes.xt[disc.mask_in_nodes] + \
                                       self.cells.left_part_size.current_value[enr_cell] / 2.
            modified_coord[index, 2] = self.nodes.xt[disc.mask_out_nodes] - \
                                       self.cells.right_part_size.current_value[enr_cell] / 2.
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
    def pressure_field(self) -> np.array:
        """
        Pressure field
        """
        return self.cells.pressure_field

    @property
    def density_field(self) -> np.array:
        """
        Density field
        """
        return self.cells.density_field

    @property
    def energy_field(self) -> np.array:
        """
        Internal energy field
        """
        return self.cells.energy_field

    @property
    def artificial_viscosity_field(self) -> np.array:
        """
        Artificial viscosity field
        """
        return self.cells.artificial_viscosity_field

    @property
    def deviatoric_stress_field(self) -> np.array:
        """
        Deviatoric stress field
        """
        return self.cells.deviatoric_stress_field

    @property
    def stress_xx_field(self) -> np.array:
        """
        First component of the Cauchy stress tensor
        """
        return self.cells.stress_xx_field
