# -*- coding: utf-8 -*-
"""
Definition of CohesiveZoneLaw class
"""
import numpy as np


class CohesiveLaw:
    """
    A class for cohesive law implementation
    """
    def __init__(self, cohesive_law_points: np.array):
        """
        Build a cohesive zone model object

        :param cohesive_law_points: array describing the stress - opening curve of the
        cohesive model
        # TODO : mettre à jour data container pour construire les modèles cohésifs
        """
        assert len(cohesive_law_points.shape) == 2, "array should be 2D"
        assert cohesive_law_points.shape[1] == 2, "array should be size (x, 2)"
        assert cohesive_law_points[0, 0] == 0., "first value of separation should be 0."
        assert cohesive_law_points[-1, 1] == 0., "last value of stress should be 0."
        self.cohesive_law_points = cohesive_law_points
        self.separation_points = self.cohesive_law_points[:, 0]
        assert np.all(np.diff(self.separation_points) >= 0), "separation is not sorted"

    def compute_cohesive_force(self, opening):
        """
        Returns the cohesive force associated with the given opening

        :param opening: discontinuity opening
        :return: float
        """
        # Theoretically, this case should not append but this verification ensure no index error
        # will occur in the future
        if opening > self.separation_points[-1]:
            return 0.

        # Find the relevant points to interpolate cohesive law
        index = np.searchsorted(self.separation_points, opening)

        # Interpolate the cohesive law
        return CohesiveLaw.interpolate_cohesive_law(opening,
                                                    self.cohesive_law_points[index - 1, 0],
                                                    self.cohesive_law_points[index, 0],
                                                    self.cohesive_law_points[index - 1, 1],
                                                    self.cohesive_law_points[index, 1])

    @classmethod
    def interpolate_cohesive_law(cls, opening, separation_1, separation_2, stress_1, stress_2):
        """
        Interpolate the value of cohesive stress between points 1 and 2

        :param opening: discontinuity opening
        :param separation_1 : separation at point 1
        :param separation_2 : separation at point 2
        :param stress_1 : stress at point 1
        :param stress_2 : stress at point 2
        :return: cohesive stress
        """
        slope = (stress_2 - stress_1) / (separation_2 - separation_1)
        return stress_1 + slope * (opening - separation_1)
