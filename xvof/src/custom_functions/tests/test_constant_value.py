# -*- coding: utf-8 -*-
"""
Unittests of the constant_pressure module
"""
import unittest

import xvof.src.boundary_condition.constant_pressure as CstPressure


class ConstantPressureTest(unittest.TestCase):
    """
    This class tests the methods of the ConstantPressure class
    """
    def setUp(self):
        """
        Tests initialization
        """
        self.value = 2.
        self.pressure = CstPressure.ConstantPressure(self.value)

    def test_evaluate(self):
        """
        Tests the evalulate method
        """
        value_test = self.pressure.evaluate(0)
        self.assertEqual(value_test, self.value)


if __name__ == '__main__':
    unittest.main()
