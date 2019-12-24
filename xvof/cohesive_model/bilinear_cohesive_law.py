# -*- coding: iso-8859-1 -*-
"""
Definition of BilinearCohesiveZoneModel
"""
from xvof.cohesive_model.cohesive_law_base import CohesiveZoneModelBase
from xvof.cohesive_model.linear_cohesive_law import LinearCohesiveZoneModel


class BilinearCohesiveZoneModel(CohesiveZoneModelBase):
    """
    An interface for all cohesive zone model
    """
    def __init__(self, cohesive_strength, separation_1, stress_1, critical_separation, unloading_model):
        """
        Constructeur :
        :param sigma_0 : cohesive strength
        :param separation_1 : ouverture au changement de pente
        :param stress_1 : contrainte au changement de pente
        :param delta_2 : ouverture critique (sigma = 0)
        """
        CohesiveZoneModelBase.__init__(self, cohesive_strength, critical_separation, unloading_model)
        self.separation_1 = separation_1
        self.stress_1 = stress_1

        # La loi bilinéaire peut être vue comme une juxtaposition de 2 lois linéaires
        self.linear_law_1 = LinearCohesiveZoneModel(stress_1=cohesive_strength, separation_1=0.,
                                                    stress_2=stress_1, separation_2=separation_1,
                                                    unloading_model=unloading_model)
        self.linear_law_2 = LinearCohesiveZoneModel(stress_1=stress_1, separation_1=separation_1,
                                                    stress_2=0., separation_2=critical_separation,
                                                    unloading_model=unloading_model)

        # Vérification de la cohérence des paramètres
        if not separation_1 <= critical_separation:
            raise ValueError("""Erreur dans le jeu de données.
            Les valeurs de séparation dans la loi cohésive ne sont pas cohérentes.
            Il faut avoir Separation Linear < Critical Separation""")

    def compute_cohesive_force_in_model(self, current_disc_opening):
        """
        Calcul de la cohesive stress dans le cas d'une loi bilinéaire
        :param current_disc_opening:
        :return:
        """
        if current_disc_opening <= self.separation_1:
            return self.linear_law_1.compute_cohesive_force_in_model(current_disc_opening)
        else:
            return self.linear_law_2.compute_cohesive_force_in_model(current_disc_opening)