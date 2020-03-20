# -*- coding: iso-8859-1 -*-
"""
Definition of LinearCohesiveZoneModel
"""

from xfv.src.cohesive_model.cohesive_law_base import CohesiveZoneModelBase


class LinearCohesiveZoneModel(CohesiveZoneModelBase):
    """
    A cohesive zone model with linear cohesive law
    """
    def __init__(self, stress_1, separation_2, unloading_model, separation_1=0., stress_2=0.):
        """
        Constructeur :
        :param stress_1 : stress point 1
        :param separation_1 : ouverture point 1
        :param stress_2 : stress point 2
        :param separation_2 : ouverture point 2
        """
        CohesiveZoneModelBase.__init__(self, stress_1, separation_2, unloading_model)
        self.stress_1 = stress_1
        self.separation_1 = separation_1
        self.stress_2 = stress_2
        self.separation_2 = separation_2

        # Vérification de la cohérence des paramètres
        if not separation_1 <= separation_2:
            raise ValueError("""Erreur dans le jeu de données. Les valeurs de séparation dans la loi cohésive ne
            sont pas cohérentes. Il faut avoir Separation 1 < Separation 2""")

    def compute_cohesive_force_in_model(self, current_disc_opening):
        """
        Calcule par une interpolation linéaire la contrainte correspondant à l'ouverture en argument
        ____
            \
             \
              _____
        :param current_disc_opening: ouverture de la discontinuité à calculer
        :return: contrainte
        """
        # if current_disc_opening <= self.separation_1:
        #     return self.stress_1
        #
        # if current_disc_opening >= self.separation_2:
        #     return self.stress_2

        return self.stress_1 + (self.stress_2 - self.stress_1) / (self.separation_2 - self.separation_1) * (current_disc_opening - self.separation_1)
