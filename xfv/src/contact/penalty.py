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

    def compute_contact_force(self, node_velocity: np.array, disc: Discontinuity,
                              delta_t: float) -> float:
        """
        Checks if contact and apply correction
        :param node_velocity: node velocity array
        :param disc: discontinuity to examine
        :param delta_t : time step
        """
        opening = disc.discontinuity_opening.new_value
        contact_force = self._compute_penalty_force(opening)
        return contact_force

    def _compute_penalty_force(self, opening) -> float:
        """
        Compute the penalty force to apply in order to penalize contact
        :param opening: discontinuity opening
        :return force = float
        """
        if opening > 0.:
            return 0.
        return self.penalty_stiffness * opening
