# -*- coding: utf-8 -*-
"""
Definition of Contact management
"""
from abc import abstractmethod
import numpy as np

from xfv.src.discontinuity.discontinuity import Discontinuity


class ContactBase(object):
    """
    An interface for all cohesive zone model
    """
    def __init__(self, *args, **kwargs):
        """
        Build the contact model to prevent discontinuities boundaries interpenetration
        """
        pass

    @abstractmethod
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
        pass
