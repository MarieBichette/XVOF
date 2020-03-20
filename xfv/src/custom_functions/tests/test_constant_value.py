# -*- coding: utf-8 -*-
"""
Unittests of the constant_pressure module
"""
import unittest

from xvof.src.custom_functions.constant_value import ConstantValue


class ConstantValueTest(unittest.TestCase):
    """
    This class tests the methods of the ConstantValue class
    """
    def setUp(self):
        """
        Tests initialization
        """
        self.value = 2.
        self.function = ConstantValue(self.value)

    def test_evaluate(self):
        """
        Tests the evalulate method
        """
        value_test = self.function.evaluate(0)
        self.assertEqual(value_test, self.value)


if __name__ == '__main__':
    unittest.main()
