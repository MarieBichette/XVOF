# -*- coding: utf-8 -*-
"""
Definition of Contact management using penalty method
"""
import numpy as np
from xfv.src.contact.contact_base import ContactBase
from xfv.src.discontinuity.discontinuity import Discontinuity


class PenaltyContact(ContactBase):
    """
    A class to manage contacts using penalty method
    """
    def __init__(self, penalty_stiffness):
        """
        Constructor
        :param penalty_stiffness: stiffness of the penalty method
        """
        super(PenaltyContact, self).__init__()
        self.penalty_stiffness = penalty_stiffness

    def compute_contact_force(self, node_velocity: np.array, disc: Discontinuity, delta_t: float):
        """
        Checks if contact and apply correction
        :param node_coord: node coordinates array
        :param node_velocity: node velocity array
        :param disc: discontinuity to examine
        :param delta_t : time step
        """
        opening = disc.discontinuity_opening.new_value
        contact_force = self._compute_penalty_force(opening)
        contact_force = self._apply_penalty_upper_bound(disc, node_velocity, contact_force, delta_t)
        return contact_force

    def _compute_penalty_force(self, opening):
        """
        Compute the penalty force to apply in order to penalize contact
        :param opening: discontinuity opening
        :return force = float
        """
        if opening > 0.:
            return 0.
        return self.penalty_stiffness * opening

    def _apply_penalty_upper_bound(self, disc: Discontinuity, node_velocity: np.array,
                                   contact_force: float, delta_t: float):
        """
        Apply an upper bound on the computed contact force in order to ensure that the
        force does not induce snapback of the contact nodes
        (nodes enter in contact and separate in the same time step)
        Computed from Non Linear Finite Element for Continua And Structure, T. Belytschko
        :param disc: discontinuity to examine
        :param node_velocity: node velocity array
        :param contact_force: force to apply
        :param delta_t : time step
        :return:
        """
        epsilon = disc.discontinuity_position

        # Compute an equivalent node mass for the discontinuity boundaries from the nodal masses
        mass_matrix_left = disc.mass_matrix_enriched.get_mass_matrix_left()
        mass_matrix_right = disc.mass_matrix_enriched.get_mass_matrix_right()
        node_mass_1 = mass_matrix_left[0, 0]
        enr_node_mass_2 = mass_matrix_left[1, 1]
        node_mass_2 = mass_matrix_right[2, 2]
        enr_node_mass_1 = mass_matrix_right[3, 3]
        mass_g = node_mass_1 * enr_node_mass_2 /\
                 (enr_node_mass_2 * (1 - epsilon) + node_mass_1 * epsilon)
        mass_d = enr_node_mass_1 * node_mass_2 / \
                 (node_mass_2 * (1 - epsilon) + enr_node_mass_1 * epsilon)

        # Compute the velocity of discontinuity boundaries
        velocity_g = (1-epsilon) * node_velocity[disc.mask_in_nodes] \
                     + epsilon * disc.additional_dof_velocity_new[0]
        velocity_d = (1-epsilon) * disc.additional_dof_velocity_new[1] \
                     + epsilon * node_velocity[disc.mask_out_nodes]

        # Compute the upper bound of penalty force
        upper_bound = mass_g * mass_d * (velocity_g - velocity_d) / (delta_t * (mass_g + mass_d))
        bounded_force = min(contact_force, upper_bound)
        return bounded_force
