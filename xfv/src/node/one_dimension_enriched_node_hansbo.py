#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
A class defining one dimension enriched nodes
"""
import numpy as np
from xfv.src.node.one_dimension_node import OneDimensionNode
from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.mass_matrix.mass_matrix_utilities import multiplication_masse


class OneDimensionHansboEnrichedNode(OneDimensionNode):
    """
    A class for the enriched nodes with Hansbo enrichment
    """

    def __init__(self, nbr_of_nodes: int, initial_positions: np.array,
                 initial_velocities: np.array, section=1.):
        """
        Build the class OneDimensionHansboEnrichedNode

        :param nbr_of_nodes: number of nodes
        :param initial_positions: initial coordinates of nodes
        :param initial_velocities: initial velocities of nodes
        :param section: section of the bar
        """
        super().__init__(nbr_of_nodes, initial_positions, initial_velocities, section=section)
        self._v_field = np.copy(self._upundemi)

    @property
    def enrichment_concerned(self):
        """
        By definition, the nodes concerned by enrichment have non null enriched shape function on
        their support

        :return: boolean mask indicating which nodes are concerned by enrichment
        """
        return ~ self.classical

    @property
    def enrichment_not_concerned(self):
        """
        :return: boolean mask indicating which nodes are not concerned by enrichment
        """
        return ~ self.enrichment_concerned

    @property
    def classical(self):
        """
        :return: boolean mask indicating which nodes are classical
        """
        return self._classical

    @property
    def enriched(self):
        """
        :return: boolean mask indicating which nodes are enriched
        """
        return ~ self.classical

    @property
    def velocity_field(self):
        """
        Accessor on the true node velocity field
        """
        return self._v_field

    def compute_complete_velocity_field(self):
        """
        Compute the true field of node velocity
        """
        self._v_field = np.copy(self._upundemi)

    def compute_enr_new_velocity(self, disc, delta_t):
        """
        Compute the new velocity enriched degree of freedom

        :param disc: the current discontinuity
        :param delta_t: float, time step
        """
        inv_matrix = disc.mass_matrix_enriched.inverse_enriched_mass_matrix_enriched_dof
        disc.enr_velocity_new[:] = \
            disc.enr_velocity_current[:] + delta_t * \
            multiplication_masse(inv_matrix, disc.enr_force)

    def coupled_enrichment_terms_compute_new_velocity(self, disc, delta_t):
        """
        Compute the coupled terms between classical and enriched dof due to non diagonal
        complete mass matrix. Takes into account nodes concerned by enrichment and
        not only the enriched nodes

        :param disc: the current discontinuity
        :param delta_t: time step
        """
        inv_matrix = disc.mass_matrix_enriched.inverse_enriched_mass_matrix_coupling_dof
        node_in = np.where(disc.mask_in_nodes)[0][0]
        node_out = np.where(disc.mask_out_nodes)[0][0]
        mask_disc = [node_in, node_out]
        disc.enr_velocity_new += np.dot(inv_matrix.transpose(),
                                        self._force[mask_disc]) * delta_t
        self._upundemi[mask_disc] += np.dot(inv_matrix, disc.enr_force) * delta_t

    def infos(self, index):
        """
        Print information
        """
        OneDimensionNode.infos(self, index)
        message = "==> classical velocity at t-1/2 = {}\n". \
            format(self.umundemi[index])
        message += "==> classical velocity at t+1/2 = {}\n". \
            format(self.upundemi[index])
        message += "==> classical force = {}\n". \
            format(self.force[index])
        message += "------------"
        for disc in [d for d in Discontinuity.discontinuity_list()
                     if index in [d.mask_in_nodes, d.mask_out_nodes]]:
            recal_ind = (disc.mask_out_nodes == index)
            message += "==> additional degree of freedom for velocity at t-1/2 = {}\n". \
                format(disc.enr_velocity_current[recal_ind])
            message += "==> additional degree of freedom for velocity at t+1/2 = {}\n". \
                format(disc.enr_velocity_current[recal_ind])
            message += "==> additional degree of freedom for force = {}\n". \
                format(disc.enr_force.current_value[recal_ind])
        print(message)

    def initialize_additional_node_dof(self, disc: Discontinuity):
        """
        Initialialise les ddl enrichis aux noeuds

        :param disc: Discontinuity
        """
        # Warning : enr_node_2 (2-) has the same velocity / coordinates as node 2 (node out)
        # Warning : enr_node_1 (1+) has the same velocity / coordinates as node 1 (node in)
        # Consequence => initialization with array is impossible => node by node initialization

        # Velocity
        disc.enr_velocity_current[0] = np.copy(self.umundemi[disc.mask_out_nodes])  # 2-
        disc.enr_velocity_current[1] = np.copy(self.umundemi[disc.mask_in_nodes])  # 1+
        disc.enr_velocity_new[0] = np.copy(self.upundemi[disc.mask_out_nodes])  # 2-
        disc.enr_velocity_new[1] = np.copy(self.upundemi[disc.mask_in_nodes])  # 1+
        # Coordinates
        disc.enr_coordinates_current[0] = np.copy(self.xt[disc.mask_out_nodes])  # 2-
        disc.enr_coordinates_current[1] = np.copy(self.xt[disc.mask_in_nodes])  # 1+
        disc.enr_coordinates_new[0] = np.copy(self.xtpdt[disc.mask_out_nodes])  # 2-
        disc.enr_coordinates_new[1] = np.copy(self.xtpdt[disc.mask_in_nodes])  # 1+

    def reinitialize_kinematics_after_contact(self, disc: Discontinuity):
        """
        Set the new velocity to the old one to cancel the increment that has lead to contact

        :param disc: discontinuity to be considered
        """
        self._upundemi[disc.mask_disc_nodes] = np.copy(self._umundemi[disc.mask_disc_nodes])
        self._xtpdt[disc.mask_disc_nodes] = np.copy(self._xt[disc.mask_disc_nodes])

    @staticmethod
    def enriched_nodes_compute_new_coordinates(disc: Discontinuity, delta_t: float):
        """
        Compute the new nodes coordinates after enrichment

        :param disc: current discontinuity
        :param delta_t: time step
        """
        disc.enr_coordinates_new[:] = disc.enr_coordinates_current[:] + delta_t * disc.enr_velocity_new[:]

    def compute_enriched_nodes_new_force(self, contrainte_xx: np.array, enr_contrainte_xx):
        """
        Compute the enriched force on enriched nodes and apply correction for classical
        force on enriched nodes (classical ddl)

        :param contrainte_xx: vecteur contrainte xx, array de taille (nb_cell, 1)
        :param enr_contrainte_xx: vecteur contrainte xx enrichie, array de taille (nb_cell, 1)
        """
        for disc in Discontinuity.discontinuity_list():
            # For each discontinuity, compute the contribution of the cracked cell to the classical
            # node forces and compute the enriched node forces for the enriched nodes of
            # the discontinuity
            cell = disc.get_ruptured_cell_id
            epsilon = disc.position_in_ruptured_element

            sigma_minus = contrainte_xx[cell]
            sigma_plus = enr_contrainte_xx[cell]

            f_node_left_minus = sigma_minus * (1 - epsilon)
            f_node_right_plus = - sigma_plus * epsilon

            f_node_right_minus = - sigma_minus * epsilon
            f_node_left_plus = sigma_plus * (1 - epsilon)

            disc.enr_force[0] = f_node_right_minus * self.section  # F2-
            disc.enr_force[1] = f_node_left_plus * self.section  # F1+
            self._force[disc.mask_in_nodes] += f_node_left_minus * self.section  # F1-
            self._force[disc.mask_out_nodes] += f_node_right_plus * self.section  # F2+

    def compute_enriched_nodes_cohesive_forces(self, cohesive_model, stress, energy):
        """
        Compute the cohesive forces for the enriched nodes

        :param cohesive_model: cohesive model
        """
        disc_list = Discontinuity.discontinuity_list()
        if not disc_list:
            return
        nb_disc = len(disc_list)

        applied_force_arr = np.ndarray((nb_disc,))
        for ind, disc in enumerate(disc_list):  
            cohesive_stress = cohesive_model.compute_cohesive_stress(disc)
            disc.cohesive_force.new_value = cohesive_stress
            disc.dissipated_energy.new_value *= self.section
            applied_force_arr[ind] = self.section * cohesive_stress

        self.apply_force_on_discontinuity_boundaries_arr(applied_force_arr)

    def apply_force_on_discontinuity_boundaries_arr(self, force: np.ndarray) -> None:
        """
        Transport the force to apply on discontinuity boundaries on the classical and enriched nodes

        :param force: value of the force to apply
        """
        epsilon_arr = Discontinuity.discontinuity_position
        f1 = ((1. - epsilon_arr).T * force).T  # F1 # pylint:disable=invalid-name
        f2 = (epsilon_arr.T * force).T  # F2 # pylint:disable=invalid-name
        Discontinuity.enr_force[:, 0] += f2  # F2-
        Discontinuity.enr_force[:, 1] -= f1  # F1+

        nodes_in = Discontinuity.in_nodes.flatten()
        nodes_out = Discontinuity.out_nodes.flatten()
        nb_disc = len(Discontinuity.discontinuity_list())
        for ind in range(nb_disc):
            self._force[nodes_in[ind]] += f1[ind]  # F1-
            self._force[nodes_out[ind]] -= f2[ind]  # F2+
