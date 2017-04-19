#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-

"""
Classe de test du module ConstantPressure
"""
import numpy as np
import unittest

import xvof.pressurelaw.constantpressure as CstPressure

class ConstantPressureTest(unittest.TestCase):
    """
    Test case utilisé pour test les fonctions du module 'ConstantPressure'
    """
    def setUp(self):
        """
        Initialisation des tests
        """
        self.value = 2.
        self.pressure = CstPressure.ConstantPressure(self.value)



    def test_evaluate(self):
        """
        Teste la fonction evaluate de ConstantPressure
        """
        value_test = self.pressure.evaluate(0)
        self.assertEqual(value_test, self.value)



if __name__ == '__main__':
    unittest.main()


