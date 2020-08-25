# -*- coding: utf-8 -*-
"""
Definition of UnloadingModelBase
"""
from abc import ABCMeta, abstractmethod


class UnloadingModelBase(metaclass=ABCMeta):  # pylint: disable=too-few-public-methods
    """
    A model for unloading reloading path in cohesive zone model
    """
    def __init__(self, *args, **kwargs):
        """
        Constructor
        """

    @abstractmethod
    def compute_unloading_reloading_condition(self, disc, new_opening):
        """
        Compute the cohesive stress in case of unloading or reloading condition
        (new_opening is less than the discontinuity maximal opening

        :param disc: discontinuity
        :param new_opening: opening of the discontinuity
        :return: cohesive stress (float)
        """
