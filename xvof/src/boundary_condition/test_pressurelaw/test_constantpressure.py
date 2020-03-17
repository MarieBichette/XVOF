#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-

"""
Classe de test du module ConstantPressure
"""
import unittest

import xvof.boundary_condition.constantpressure as CstPressure


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
        print __name__ + " : Test evaluate value"
        value_test = self.pressure.evaluate(0)
        self.assertEqual(value_test, self.value)
        print "__[OK]"


if __name__ == '__main__':
    unittest.main()


