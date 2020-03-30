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
                      disc: Discontinuity, delta_t: float):
        """
        Checks if contact and apply correction
        :param node_coord: node coordinates array
        :param node_velocity: node velocity array
        :param node_force : node force array
        :param disc: discontinuity to examine
        :param delta_t : time step
        """
        is_contact = self._check_contact(node_coord, disc)
        if is_contact:
            contact_force = self._compute_penalty_force()
            self._apply_contact_force(disc, node_force, contact_force)

    def _check_contact(self, node_coord, disc):
        """

        :param node_coord: node coordinates array
        :param disc: discontinuity to examine
        :return:
        """
        return is_contact

    def _compute_penalty_force(self):
        """
        """
        return contact_force

    def _apply_contact_force(self, disc: Discontinuity, node_force: np.array, contact_force: float):
        """
        :param disc: discontinuity to examine
        :param node_force : node force array
        :param contact_force: force to apply
        :return:
        """
