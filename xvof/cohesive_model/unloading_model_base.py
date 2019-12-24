# -*- coding: iso-8859-1 -*-
"""
Definition of UnloadingModelBase
"""
from abc import abstractmethod


class UnloadingModelBase(object):
    """
    A model for unloading reloading path in cohesive zone model
    """
    def __init__(self):
        """
        Constructor :
        """
        pass

    @abstractmethod
    def compute_unloading_reloading_condition(self, disc, new_opening):
        """
        Calcule
        :return: contrainte
        """
        pass

    @abstractmethod
    def apply_penalty_condition(self, disc, new_opening, ancien_opening, ancien_force):
        """

        :param disc:
        :return:
        """
        pass

