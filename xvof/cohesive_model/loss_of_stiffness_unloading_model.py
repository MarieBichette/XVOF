# -*- coding: iso-8859-1 -*-
"""
Definition of UnloadingModel with opposite force
"""
from abc import abstractmethod
from xvof.cohesive_model.unloading_model_base import UnloadingModelBase


class LossOfStiffnessUnloadingModel(UnloadingModelBase):
    """
    A model for unloading reloading path in cohesive zone model
    """
    def __init__(self, slope, cohesive_strength):
        """
        Constructor :
        """
        super(LossOfStiffnessUnloadingModel, self).__init__()
        self.slope = slope
        self.sigma_0 = cohesive_strength

    @abstractmethod
    def compute_unloading_reloading_condition(self, disc, new_opening):
        """
        Calcule
        :return: contrainte
        """
        # Correction :
        cohesive_force = 0

        if disc.history_max_opening > 1.e-16:
            cohesive_force = disc.history_min_cohesive_force + disc.history_min_cohesive_force / disc.history_max_opening * (new_opening - disc.history_max_opening)

        return cohesive_force

    @abstractmethod
    def apply_penalty_condition(self, disc, new_opening, ancien_opening, ancien_force):
        """

        :param disc:
        :return:
        """
        # cohesive_force = ancien_force + 20 * ancien_force * (new_opening / ancien_opening - 1)
        cohesive_force = ancien_force + 1.e+10 * (new_opening - ancien_opening)
        return cohesive_force

