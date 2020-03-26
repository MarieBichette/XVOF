# -*- coding: iso-8859-1 -*-
"""
Definition of UnloadingModel with opposite force
"""
from abc import abstractmethod
from xfv.src.cohesive_model.unloading_model_base import UnloadingModelBase


class ZeroForceUnloadingModel(UnloadingModelBase):
    """
    A model for unloading reloading path in cohesive zone model
    """
    def __init__(self, slope, cohesive_strength):
        """
        Constructor :
        """
        super(ZeroForceUnloadingModel, self).__init__()
        self.sigma_0 = cohesive_strength
        self.slope = slope

    @abstractmethod
    def compute_unloading_reloading_condition(self, disc, new_opening):
        """
        Calcule
        :return: contrainte
        """
        return 0.

    @abstractmethod
    def apply_penalty_condition(self, disc, new_opening, ancien_opening, ancien_force):
        """

        :param disc:
        :return:
        """
        cohesive_force = ancien_force + 20 * ancien_force / ancien_opening * (new_opening - ancien_opening)
        return cohesive_force

