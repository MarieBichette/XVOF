#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-

"""
Classe de test du module RampPressure
"""
import unittest

import xvof.src.boundary_condition.creneau_ramp_pressure as CrRpPressure

class CreneauRampPressureTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'RampPressure'
    """
    def setUp(self):
        """
        Initialisation des tests
        """
        self._first_value = 2.
        self._second_value = 4.
        self._third_value = 0.
        self._start_time = 1.
        self._end_first_slope_time = 2.
        self._end_creneau = 3.
        self._end_time = 4.
        self._pressure = CrRpPressure.CreneauRampPressure(self._first_value, self._second_value, self._third_value,
                                                          self._start_time, self._end_first_slope_time,
                                                          self._end_creneau, self._end_time)

    def test_evaluate_first_value(self):
        """
        Teste la fonction evaluate de CreneauRampPressure pour t < start time
        """
        print __name__ + " : Test evaluate_first_value"
        value_test = self._pressure.evaluate(self._start_time - 10.)
        self.assertEqual(value_test, self._first_value)
        print "__[OK]"

    def test_evaluate_third_value(self):
        """
        Teste la fonction evaluate de CreneauRampPressure  pour t > end time
        """
        print __name__ + " : Test evaluate third value"
        value_test = self._pressure.evaluate(self._end_time + 10.)
        self.assertEqual(value_test, self._third_value)
        print "__[OK]"

    def test_evaluate_medium_value_first(self):
        """
        Teste la fonction evaluate de CreneauRampPressure  pour t dans la première pente du créneau
        """
        print __name__ + " : Test evaluate medium value"
        value_test = self._pressure.evaluate(0.5 * (self._start_time + self._end_first_slope_time))
        self.assertEqual(value_test, 0.5 * (self._first_value + self._second_value))
        print "__[OK]"

    def test_evaluate_creneau_value(self):
        """
        Teste la fonction evaluate de CreneauRampPressure  pour t dans la partie constante du créneau
        """
        print __name__ + " : Test creneau value"
        value_test = self._pressure.evaluate((self._end_first_slope_time + self._end_creneau) * 0.5)
        self.assertEqual(value_test, self._second_value)
        print "__[OK]"

    def test_evaluate_medium_value_second(self):
        """
        Teste la fonction evaluate de CreneauRampPressure  pour t dans la deuxième pente du créneau
        """
        print __name__ + " : Test medium second value"
        value_test = self._pressure.evaluate(0.5 * (self._end_creneau + self._end_time))
        self.assertEqual(value_test, 0.5 * (self._second_value + self._third_value))
        print "__[OK]"


if __name__ == '__main__':
    unittest.main()