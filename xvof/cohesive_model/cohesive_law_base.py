# -*- coding: iso-8859-1 -*-
"""
Definition of CohesiveZoneModelBase interface
"""
from abc import abstractmethod


class CohesiveZoneModelBase(object):
    """
    An interface for all cohesive zone models
    """
    def __init__(self, cohesive_strength, critical_separation, unloading_model):
        """
        Construction d'un modèle cohésif
        :param cohesive_strength:
        :param critical_separation:
        """
        self.cohesive_strength = cohesive_strength
        self.critical_separation = critical_separation
        self.unloading_model = unloading_model

    def compute_cohesive_stress(self, disc):
        """
        Compute the cohesive force for the current opening of discontinuity (selon les cas)
        """
        cohesive_force = 0.
        new_opening = disc.discontinuity_opening.new_value[0]

        # if new_opening < 0:
        #     # on passe plus dans ce service car on a ajouté un package pour gérer le contact
        #     print "Apply penalty condition to avoid overlapping for discontinuity {:}".format(disc.label)
        #     cohesive_force = self.apply_penalty_condition(disc, new_opening)
        #     return cohesive_force
        #     # return 0

        if disc.damage_variable.current_value[0] < 1:

            if new_opening < disc.history_max_opening:
                cohesive_force = self.compute_unloading_reloading_condition(disc, new_opening)

            elif new_opening >= disc.history_max_opening and new_opening < self.critical_separation:
                cohesive_force = self.compute_cohesive_force_in_model(new_opening)
                disc.history_max_opening = max(disc.history_max_opening, new_opening)
                disc.history_min_cohesive_force = self.compute_cohesive_force_in_model(disc.history_max_opening)

            if new_opening >= self.critical_separation:
                # print "Discontinuity " + str(disc.label) + "is completely open"
                disc.damage_variable.new_value = 1.
                cohesive_force = 0.
                disc.history_max_opening = max(disc.history_max_opening, new_opening)
                disc.history_min_cohesive_force = min(disc.history_min_cohesive_force, cohesive_force)

        return cohesive_force

    @abstractmethod
    def compute_cohesive_force_in_model(self, opening):
        """
        Le calcul est délégué à la loi cohésive choisie
        :param opening:
        :return:
        """
        return 0.

    def apply_penalty_condition(self,  disc, new_opening):
        """
        Condition de pénalité pour éviter que l'ouverture devienne négative
        :param disc:
        :return:
        """
        # ancien_opening = disc.discontinuity_opening.current_value[0]
        # ancien_force = disc.cohesive_force.current_value[0]
        ancien_opening = 0
        ancien_force = self.cohesive_strength
        return self.unloading_model.apply_penalty_condition(disc, new_opening, ancien_opening,
                                                                          ancien_force)
        # return 0.

    def compute_unloading_reloading_condition(self, disc, new_opening):
        """
        Charge / Décharge de la zone cohésive quand on est en dessous de l'ouverture max atteinte
        :param opening: ouverture courante
        :param disc:
        :return:
        """
        # ancien_opening = disc.discontinuity_opening.current_value[0]
        # ancien_force = disc.cohesive_force.current_value[0]
        return self.unloading_model.compute_unloading_reloading_condition(disc, new_opening)