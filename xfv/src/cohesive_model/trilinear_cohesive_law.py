# -*- coding: iso-8859-1 -*-
"""
Definition of TrilinearCohesiveZoneModel
"""
from xfv.src.cohesive_model.cohesive_law_base import CohesiveZoneModelBase
from xfv.src.cohesive_model.linear_cohesive_law import LinearCohesiveZoneModel


class TrilinearCohesiveZoneModel(CohesiveZoneModelBase):
    """
    An interface for all cohesive zone model
    """
    def __init__(self, cohesive_strength, separation_1, stress_1, separation_2, stress_2, critical_separation,
                 unloading_model):
        """
        Constructeur :
        :param cohesive_strength : cohesive strength
        :param separation_1 : ouverture au changement de pente 1
        :param stress_1 : contrainte au changement de pente 1
        :param separation_2 : ouverture au changement de pente 2
        :param stress_2 : contrainte au changement de pente 2
        :param critical_separation : ouverture critique (sigma = 0)
        """
        CohesiveZoneModelBase.__init__(self, cohesive_strength, critical_separation, unloading_model)
        self.separation_1 = separation_1
        self.stress_1 = stress_1
        self.separation_2 = separation_2
        self.stress_2 = stress_2

        # La loi trilinéaire peut être vue comme une juxtaposition de 3 lois linéaires
        self.linear_law_1 = LinearCohesiveZoneModel(stress_1=cohesive_strength, separation_1=0.,
                                                    stress_2=stress_1, separation_2=separation_1,
                                                    unloading_model=unloading_model)
        self.linear_law_2 = LinearCohesiveZoneModel(stress_1=stress_1, separation_1=separation_1,
                                                    stress_2=stress_2, separation_2=separation_2,
                                                    unloading_model=unloading_model)
        self.linear_law_3 = LinearCohesiveZoneModel(stress_1=stress_2, separation_1=separation_2,
                                                    stress_2=0., separation_2=critical_separation,
                                                    unloading_model=unloading_model)

        # Vérification de la cohérence des paramètres
        if separation_1 > separation_2 or separation_2 > critical_separation:
            raise ValueError("""Erreur dans le jeu de données.
            Les valeurs de séparation dans la loi cohésive ne sont pas cohérentes.
            Il faut avoir Separation 1 < Separation 2 < Critical Separation""")

    def compute_cohesive_force_in_model(self, current_disc_opening):
        """
        Calcul de la cohesive stress dans le cas d'une loi trilinéaire
        :param current_disc_opening:
        :return:
        """
        if current_disc_opening <= self.separation_1:
            return self.linear_law_1.compute_cohesive_force_in_model(current_disc_opening)
        elif current_disc_opening <= self.separation_2:
            return self.linear_law_2.compute_cohesive_force_in_model(current_disc_opening)
        else:
            return self.linear_law_3.compute_cohesive_force_in_model(current_disc_opening)
