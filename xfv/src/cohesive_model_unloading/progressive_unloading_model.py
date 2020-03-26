# -*- coding: iso-8859-1 -*-
"""
Definition of UnloadingModel with opposite force
"""
from abc import abstractmethod
from xfv.src.cohesive_model_unloading.unloading_model_base import UnloadingModelBase


class ProgressiveUnloadingModel(UnloadingModelBase):
    """
    A model for unloading reloading path in cohesive zone model
    """
    def __init__(self, slope, cohesive_strength):
        """
        Constructor :
        """
        super(ProgressiveUnloadingModel, self).__init__()
        self.slope = slope
        self.sigma_0 = cohesive_strength

    @abstractmethod
    def compute_unloading_reloading_condition(self, disc, new_opening):
        """
        Calcule
        :return: contrainte
        """
        cohesive_force = 0
        if disc.history_max_opening > 1.e-15:
            cohesive_force = disc.history_min_cohesive_force + self.slope * (new_opening - disc.history_max_opening)
        return cohesive_force

    @abstractmethod
    def apply_penalty_condition(self, disc, new_opening, ancien_opening, ancien_force):
        """

        :param disc:
        :return:
        """
        cohesive_force = ancien_force + 1000 * self.slope * (new_opening - ancien_opening)
        return cohesive_force

