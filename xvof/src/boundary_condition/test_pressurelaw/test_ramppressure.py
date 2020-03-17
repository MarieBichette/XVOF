#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-

"""
Classe de test du module RampPressure
"""
import unittest
import xvof.src.boundary_condition.ramppressure as RpPressure


class RampPressureTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'RampPressure'
    """
    def setUp(self):
        """
        Initialisation des tests
        """
        self._first_value = 2.
        self._second_value = 4.
        self._start_time = 1.
        self._end_time = 3.
        self._pressure = RpPressure.RampPressure(self._first_value, self._second_value, self._start_time, self._end_time)


    def test_evaluate_first_value(self):
        """
        Teste la fonction evaluate de RampPressure pour t < start time
        """
        print __name__ + " : Test evaluate first value"
        value_test = self._pressure.evaluate(self._start_time - 10.)
        self.assertEqual(value_test, self._first_value)
        print "__[OK]"
        

    def test_evaluate_second_value(self):
        """
        Teste la fonction evaluate de RampPressure  pour t > end time
        """
        print __name__ + " : Test evaluate second value"
        value_test = self._pressure.evaluate(self._end_time + 10.)
        self.assertEqual(value_test, self._second_value)
        print "__[OK]"

    def test_evaluate_medium_value(self):
        """
        Teste la fonction evaluate de RampPressure  pour t > start time + end time / 2
        """
        print __name__ + " : Test evaluate medium value"
        value_test = self._pressure.evaluate(0.5 * (self._start_time + self._end_time))
        self.assertEqual(value_test, 0.5 * (self._first_value + self._second_value))
        print "__[OK]"


if __name__ == '__main__':
    unittest.main()