# -*- coding: utf-8 -*-
"""
Definition of Contact management using penalty method
"""
from xfv.src.contact.contact_base import ContactBase
from xfv.src.discontinuity.discontinuity import Discontinuity


class PenaltyContact(ContactBase):
    """
    An interface for all cohesive zone model
    """
    def __init__(self, penalty_stiffness):
        """
        Constructor
        :param penalty_stiffness: stiffness of the penalty method
        """
        super(PenaltyContact, self).__init__()
        self.penalty_stiffness = penalty_stiffness

    def apply_contact(self, node_coord: np.array, node_velocity: np.array, node_force: np.array,
                      disc: Discontinuity, delta_t: float, section: float):
        """
        Checks if contact and apply correction
        :param node_coord: node coordinates array
        :param node_velocity: node velocity array
        :param node_force : node force array
        :param disc: discontinuity to examine
        :param delta_t : time step
        :param section : section of the bar
        """
        opening = self._compute_discontinuity_opening(node_coord, disc)
        contact_force = self._compute_penalty_force(opening)
        contact_force = self._apply_penalty_upper_bound(contact_force, delta_t)
        self._apply_contact_force(disc, node_force, contact_force, section)

    @staticmethod
    def _compute_discontinuity_opening(self, node_coord, disc):
        """
        Compute the discontinuity opening
        :param node_coord: node coordinates array
        :param disc: discontinuity to examine
        :return: discontinuity opening (float)
        """
        left_boundary = node_coord[disc.mask_in_nodes] + disc.left_size.new_value
        right_boundary = node_coord[disc.mask_out_nodes] - disc.right_size.new_value
        return right_boundary - left_boundary

    def _compute_penalty_force(self, opening):
        """
        Compute the penalty force to apply in order to penalize contact
        :param opening: discontinuity opening
        :return force = float
        """
        if opening > 0.:
            return 0.
        return self.penalty_stiffness * opening

    def _apply_penalty_upper_bound(self, contact_force: float, delta_t: float):
        """
        Apply an upper bound on the computed contact force in order to ensure that the
        force does not induce snapback of the contact nodes
        (nodes enter in contact and separate in the same time step)
        Computed from Non Linear Finite Element for Continua And Structure, T. Belytschko
        :param contact_force: force to apply
        :param delta_t : time step
        :return:
        """
        epsilon = self.disc.discontinuity_position

        # Compute an equivalent node mass for the discontinuity boundaries from the nodal masses
        mass_matrix_left = self.disc.mass_matrix_enriched.get_mass_matrix_left()
        mass_matrix_right = self.disc.mass_matrix_enriched.get_mass_matrix_right()
        node_mass_1 = mass_matrix_left[0, 0]
        enr_node_mass_2 = mass_matrix_left[1, 1]
        node_mass_2 = mass_matrix_right[2, 2]
        enr_node_mass_1 = mass_matrix_right[3, 3]
        mass_g = node_mass_1 * enr_node_mass_2 /\
                 (enr_node_mass_2 * (1 - epsilon) + node_mass_1 * epsilon)
        mass_d = enr_node_mass_1 * node_mass_2 / \
                 (node_mass_2 * (1 - epsilon) + enr_node_mass_1 * epsilon)

        # Compute the velocity of discontinuity boundaries
        velocity_g = (1-epsilon) * node_velocity[self.disc.mask_in_nodes] \
                     + epsilon * self.disc.additional_dof_velocity_new[0]
        velocity_d = (1-epsilon) * self.disc.additional_dof_velocity_new[1] \
                     + epsilon * node_velocity[self.disc.mask_out_nodes]

        # Compute the upper bound of penalty force
        upper_bound = mass_g * mass_d * (velocity_g - velocity_d) / (delta_t * (mass_g + mass_d))
        bounded_force = min(contact_force, upper_bound)
        return bounded_force

    @staticmethod
    def _apply_contact_force(self, disc: Discontinuity, node_force: np.array,
                             contact_force: float, section: float):
        """
        Apply the contact forces on classical nodal forces and enriched nodal forces
        :param disc: discontinuity to examine
        :param node_force : node force array
        :param contact_force: force to apply
        :param section: section of the bar
        :return:
        """
        f_contact = section * contact_force
        epsilon = disc.position_in_ruptured_element

        # Apply cohesive stress on enriched nodes
        node_force[disc.mask_in_nodes] += (1. - epsilon) * f_contact  # F1-
        disc.additional_dof_force[0] += epsilon * f_contact  # F2-
        node_force[disc.mask_out_nodes] += - epsilon * f_contact  # F2+
        disc.additional_dof_force[1] += - (1. - epsilon) * f_contact  # F1+
