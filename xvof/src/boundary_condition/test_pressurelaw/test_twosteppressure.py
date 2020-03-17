#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-

"""
Classe de test du module TwoStepsPressure
"""
import unittest

import xvof.src.boundary_condition.twostepspressure as TwoStpPressure

class TwoStepsPressureTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'ConstantPressure'
    """
    def setUp(self):
        """
        Initialisation des tests
        """
        self._first_value = 2.
        self._second_value = 4.
        self._critical_time = 0.
        self._pressure = TwoStpPressure.TwoStepsPressure(self._first_value, self._second_value, self._critical_time)


    def test_evaluate_fisrt_value(self):
        """
        Teste la fonction evaluate de TwoStepPressure pour t < critical time
        """
        print __name__ + " : Test evaluate first value"
        value_test = self._pressure.evaluate(self._critical_time - 10.)
        self.assertEqual(value_test, self._first_value)
        print "__[OK]"
        

    def test_evaluate_second_value(self):
        """
        Teste la fonction evaluate de TwoStepPressure pour t > critical time
        """
        print __name__ + " : Test evaluate second value"
        value_test = self._pressure.evaluate(self._critical_time + 10.)
        self.assertEqual(value_test, self._second_value)
        print "__[OK]"


if __name__ == '__main__':
    unittest.main()