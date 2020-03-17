#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-

"""
Classe de test du module ConstantPressure
"""
import numpy as np
import unittest

import xvof.pressurelaw.twostepspressure as TwoStpPressure

class ConstantPressureTest(unittest.TestCase):
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
        Teste la fonction evaluate de ConstantPressure
        """
        value_test = self._pressure.evaluate(self._critical_time - 10.)
        self.assertEqual(value_test, self._first_value)
        

    def test_evaluate_second_value(self):
        """
        Teste la fonction evaluate de ConstantPressure
        """
        value_test = self._pressure.evaluate(self._critical_time + 10.)
        self.assertEqual(value_test, self._second_value)



if __name__ == '__main__':
    unittest.main()