# -*- coding: utf-8 -*-
"""
Definition of CohesiveZoneModel class to manage the cohesive discontinuities
"""
import numpy as np

from xfv.src.cohesive_model.cohesive_law import CohesiveLaw
from xfv.src.cohesive_model_unloading.unloading_model_base import UnloadingModelBase


class CohesiveZoneModel:
    """
    A class for the computation of the cohesive force
    """
    def __init__(self, cohesive_law_points: np.array, unloading_model: UnloadingModelBase):
        """
        Construction d'un modèle cohésif
        :param cohesive_law_points: array describing the stress - opening curve of the
        cohesive model
        # TODO : mettre à jour data container pour construire les modèles cohésifs
        :param unloading_model:
        """
        self._critical_separation = cohesive_law_points[-1, 0]
        self._cohesive_law = CohesiveLaw(cohesive_law_points)
        self._unloading_model = unloading_model

    def compute_cohesive_stress(self, disc):
        """
        Compute the cohesive force for the current opening of discontinuity according to the
        current discontinuity opening
        :param disc :discontinuity
        """
        cohesive_force = 0.
        new_opening = disc.discontinuity_opening.new_value[0]

        if disc.damage_variable.current_value[0] < 1:
            if new_opening < disc.history_max_opening:
                cohesive_force = \
                    self._unloading_model.compute_unloading_reloading_condition(disc, new_opening)

            elif disc.history_max_opening <= new_opening < self._critical_separation:
                cohesive_force = self._cohesive_law.compute_cohesive_force(new_opening)
                # Update the discontinuity indicators
                disc.history_max_opening = max(abs(disc.history_max_opening), max(new_opening))
                disc.history_min_cohesive_force = \
                    self._cohesive_law.compute_cohesive_force(disc.history_max_opening)
                disc.damage_variable.new_value = new_opening / self._critical_separation

            if new_opening >= self._critical_separation:
                disc.damage_variable.new_value = 1.
                cohesive_force = 0.
                disc.history_max_opening = max(abs(disc.history_max_opening), abs(new_opening))
                disc.history_min_cohesive_force = 0.

        return cohesive_force
