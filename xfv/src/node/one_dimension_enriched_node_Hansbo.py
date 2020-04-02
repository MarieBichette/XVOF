#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
A class defining one dimension enriched nodes
"""
import numpy as np
from xfv.src.node.one_dimension_enriched_node import OneDimensionEnrichedNode
from xfv.src.node.one_dimension_node import OneDimensionNode
from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.data.data_container import DataContainer


class OneDimensionHansboEnrichedNode(OneDimensionEnrichedNode):
    """
    A class for the enriched nodes with Hansbo enrichment
    """

    def __init__(self, nbr_of_nodes, initial_positions, initial_velocities, section=1.):
        super(OneDimensionHansboEnrichedNode, self).__init__(nbr_of_nodes, initial_positions,
                                                             initial_velocities, section=section)
        self._v_field = np.copy(self._upundemi)

    @property
    def enrichment_concerned(self):
        """
        :return: boolean mask indicating which nodes are concerned by enrichment
        (fonctions de formes enrichies N* non nulles sur leur support)
        en 1D : correspond aux noeuds enrichis et les premiers voisins.
        """
        return ~ self.classical

    @property
    def enrichment_not_concerned(self):
        """
        :return: boolean mask indicating which nodes are not concerned by enrichment
        """
        return ~ self.enrichment_concerned

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
                format(disc.additional_dof_velocity_current[recal_ind])
            message += "==> additional degree of freedom for velocity at t+1/2 = {}\n". \
                format(disc.additional_dof_velocity_current[recal_ind])
            message += "==> additional degree of freedom for force = {}\n". \
                format(disc.additional_dof_force.current_value[recal_ind])
        print(message)

    def initialize_additional_node_dof(self, disc: Discontinuity):
        """
        Initialialise les ddl enrichis aux noeuds
        :param disc: Discontinuity
        """
        # Velocity
        disc.additional_dof_velocity_current[:] = np.copy(self.umundemi[disc.mask_disc_nodes])
        disc.additional_dof_velocity_new[:] = np.copy(self.upundemi[disc.mask_disc_nodes])
        # Coordinates
        disc.additional_dof_coordinates_current[:] = np.copy(self.xtpdt[disc.mask_disc_nodes])
        disc.additional_dof_coordinates_new[:] = np.copy(self.xt[disc.mask_disc_nodes])

    def reinitialize_kinematics_after_contact(self, disc: Discontinuity):
        """
        Set the new velocity to the old one to cancel the increment that has lead to contact
        :param disc : discontinuity to be considered
        """
        self._upundemi[disc.mask_disc_nodes] = np.copy(self._umundemi[disc.mask_disc_nodes])
        self._xtpdt[disc.mask_disc_nodes] = np.copy(self._xt[disc.mask_disc_nodes])

    def enriched_nodes_compute_new_coordinates(self, delta_t: float):
        """
        Compute the new nodes coordinates after enrichment
        :param delta_t: time step
        """
        for disc in Discontinuity.discontinuity_list():
            disc._additional_dof_coordinates_new = disc.additional_dof_coordinates_current + \
                                                   delta_t * disc.additional_dof_velocity_new

    def compute_enriched_nodes_new_force(self, contrainte_xx: np.array):
        """
        Compute the enriched force on enriched nodes and apply correction for classical
        force on enriched nodes (classical ddl)
        :param contrainte_xx : vecteur contrainte xx, array de taille (nb_cell, 1)
        """
        for disc in Discontinuity.discontinuity_list():
            # For each discontinuity, compute the contribution of the cracked cell to the classical
            # node forces and compute the enriched node forces for the enriched nodes of
            # the discontinuity
            cell = disc.ruptured_cell_id
            epsilon = disc.position_in_ruptured_element

            sigma_minus = contrainte_xx[cell]
            sigma_plus = disc.additional_dof_stress[:, 0]

            f_node_left_minus = sigma_minus * (1 - epsilon)
            f_node_right_plus = - sigma_plus * epsilon

            f_node_right_minus = - sigma_minus * epsilon
            f_node_left_plus = sigma_plus * (1 - epsilon)

            disc.additional_dof_force[0] = f_node_right_minus * self.section  # F2-
            disc.additional_dof_force[1] = f_node_left_plus * self.section  # F1+
            self._force[disc.mask_in_nodes] += f_node_left_minus * self.section  # F1-
            self._force[disc.mask_out_nodes] += f_node_right_plus * self.section  # F2+

    def compute_enriched_nodes_cohesive_forces(self, cohesive_model):
        """
        Compute the cohesive forces for the enriched nodes
        :param cohesive_model : cohesive model
        """
        for disc in Discontinuity.discontinuity_list():
            # Compute cohesive stress
            cohesive_stress = cohesive_model.compute_cohesive_stress(disc)
            disc.cohesive_force.new_value = cohesive_stress
            self.apply_force_on_discontinuity_boundaries(disc, cohesive_stress)

    def apply_force_on_discontinuity_boundaries(self, disc: Discontinuity, stress: float):
        """
        Transport the force to apply on discontinuity boundaries on the classical and enriched nodes
        :param disc: current discontinuity
        :param stress: value of the force to apply
        :return:
        """
        applied_force = self.section * stress
        epsilon = disc.position_in_ruptured_element

        # Apply cohesive stress on enriched nodes
        self._force[disc.mask_in_nodes] += (1. - epsilon) * applied_force  # F1-
        disc.additional_dof_force[0] += epsilon * applied_force  # F2-
        self._force[disc.mask_out_nodes] += - epsilon * applied_force  # F2+
        disc.additional_dof_force[1] += - (1. - epsilon) * applied_force  # F1+
