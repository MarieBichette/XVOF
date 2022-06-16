# -*- coding: utf-8 -*-
# pylint: disable=protected-access
"""
Definition of Contact management
"""
import numpy as np
from xfv.src.contact.contact_base import ContactBase
from xfv.src.discontinuity.discontinuity import Discontinuity


class LagrangianMultiplierContact(ContactBase):  # pylint: disable=too-few-public-methods
    """
    A class for contact management using Lagrangian multipliers
    """
    def __init__(self):
        """
        Constructor
        """
        super().__init__()
        self._lambda_multiplier = 0

    def compute_contact_force(self, node_velocity: np.array, disc: Discontinuity,
                              delta_t: float) -> float:
        """
        Checks if contact and apply correction

        :param node_velocity: node velocity array
        :param disc: discontinuity to examine
        :param delta_t: time step
        """
        if disc.discontinuity_opening.new_value >= 0.:
            return 0.
        # else
        return LagrangianMultiplierContact._compute_lagrangian_multiplier(node_velocity,
                                                                          disc, delta_t)

    @staticmethod
    def _compute_lagrangian_multiplier(node_velocity: np.array, disc: Discontinuity,
                                       time_step: float) -> float:
        """
        Compute the lagrangian multiplier representing the contact force to apply to ensure the
        non penetration of the discontinuity boundaries
        Computed from Non Linear Finite Element for Continua And Structure, T. Belytschko, p.610

        :param node_velocity: array of the predicted velocities
        :param disc: current discontinuity
        :param time_step: time step
        :return:
        """
        epsilon = disc.discontinuity_position
        # Computes fictive node mass for the discontinuity boundaries
        m_left = disc.mass_matrix_enriched.get_mass_matrix_left()
        m_right = disc.mass_matrix_enriched.get_mass_matrix_right()
        m_1 = m_left[0, 0]
        m_2_enr = m_left[1, 1]
        m_2 = m_right[2, 2]
        m_1_enr = m_right[3, 3]
        node_mass_g = m_1 * m_2_enr / (m_2_enr * (1 - epsilon) + m_1 * epsilon)
        node_mass_d = m_1_enr * m_2 / (m_2 * (1 - epsilon) + m_1_enr * epsilon)
        # Compute the velocities of the discontinuity boundaries
        velocity_g = ((1 - epsilon) * node_velocity[disc.in_nodes]
                      + epsilon * disc.enr_velocity_new[0])
        velocity_d = ((1 - epsilon) * disc.enr_velocity_new[1]
                      + epsilon * node_velocity[disc.out_nodes])
        # Compute the contact force in order to have ug = ud afterward
        lambda_multiplier = (velocity_d - velocity_g).flatten() / (
            time_step * (node_mass_g + node_mass_d) / (node_mass_g * node_mass_d))
        # signe inverse par rapport au papier car les forces ont un sens inverse
        return lambda_multiplier
