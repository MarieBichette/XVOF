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
    def __init__(self, cohesive_law_points: np.array, unloading_model: UnloadingModelBase, cohesive_zone_model_name, dissipated_energy, purcentage):
        """
        Construction d'un modèle cohésif

        :param cohesive_law_points: array describing the stress - opening curve of the
        cohesive model
        :param unloading_model:
        """
        self._critical_separation = cohesive_law_points[-1, 0]
        self._critical_strength = cohesive_law_points[0,-1]
        self._cohesive_law = CohesiveLaw(cohesive_law_points)
        self._unloading_model = unloading_model
        self._cohesive_zone_model_name = cohesive_zone_model_name
        self._dissipated_energy = dissipated_energy
        self._purcentage = purcentage

    @property
    def cohesive_zone_model_name(self):  # pylint: disable=invalid-name
        """
        Cohesive zone model name
        """
        return self._cohesive_zone_model_name

    @property
    def critical_separation(self):  # pylint: disable=invalid-name
        """
        critical separation
        """
        return self._critical_separation

    @property
    def critical_strength(self):  # pylint: disable=invalid-name
        """
        critical strength
        """
        return self._critical_strength

    @property
    def dissipated_energy(self):  # pylint: disable=invalid-name
        """
        dissipated energy
        """
        return self._dissipated_energy

    @property
    def purcentage(self):  # pylint: disable=invalid-name
        """
        percentage of internal energy to be dissipated
        """
        return self._purcentage

    def compute_cohesive_stress(self, disc):
        """
        Compute the cohesive force for the current opening of discontinuity according to the
        current discontinuity opening

        :param disc: discontinuity

	    Attention : formule de la dissipation d'energie du modèle cohesif seulement dans le cas de la décharge à zéro
        """
        cohesive_force = 0.
        new_opening = disc.discontinuity_opening.new_value[0]
        self._cohesive_law = CohesiveLaw(np.array([[0, disc.critical_strength], [disc.critical_separation, 0]]))

        if new_opening <= 0. and disc.history_max_opening < disc.critical_separation :
            cohesive_force = 0.
            disc.dissipated_energy.new_value = disc.critical_strength*disc.history_max_opening/2.

        elif 0. < new_opening < disc.history_max_opening and disc.history_max_opening < disc.critical_separation :
            cohesive_force = \
                self._unloading_model.compute_unloading_reloading_condition(disc, new_opening)
            disc.dissipated_energy.new_value = disc.critical_strength*disc.history_max_opening/2.

        elif disc.history_max_opening <= new_opening < disc.critical_separation :
            cohesive_force = self._cohesive_law.compute_cohesive_force(new_opening)
            # Update the discontinuity indicators
            disc.history_max_opening = max(abs(disc.history_max_opening), abs(new_opening))
            disc.history_min_cohesive_force = \
                self._cohesive_law.compute_cohesive_force(disc.history_max_opening)
            disc.damage_variable.new_value = new_opening / disc.critical_separation
            disc.dissipated_energy.new_value = disc.critical_strength*disc.history_max_opening/2.
        else :
            disc.damage_variable.new_value = 1.
            cohesive_force = 0.
            disc.history_max_opening = max(abs(disc.history_max_opening), abs(new_opening))
            disc.history_min_cohesive_force = 0.
            disc.dissipated_energy.new_value = disc.critical_strength*disc.critical_separation/2.

        return cohesive_force
