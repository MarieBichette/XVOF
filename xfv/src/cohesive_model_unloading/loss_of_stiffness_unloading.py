#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Definition of UnloadingModel with loss of stiffness
"""
from xfv.src.cohesive_model_unloading.unloading_model_base import UnloadingModelBase


class LossOfStiffnessUnloading(UnloadingModelBase):
    """
    A model for unloading reloading path with decreasing stiffness
    """
    def __init__(self):
        """
        Constructor
        """
        super(LossOfStiffnessUnloading, self).__init__()

    def compute_unloading_reloading_condition(self, disc, new_opening):
        """
        Compute the cohesive stress in case of unloading or reloading condition
        (new_opening is less than the discontinuity maximal opening
        :param disc : discontinuity
        :param new_opening : opening of the discontinuity
        :return: cohesive stress (float)
        """
        slope = disc.history_min_cohesive_force / disc.history_max_opening
        cohesive_force = disc.history_min_cohesive_force +\
                         slope * (new_opening - disc.history_max_opening)
        return cohesive_force
